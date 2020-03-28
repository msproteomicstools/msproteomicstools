#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2013, ETH Zurich
For a full list of authors, refer to the file AUTHORS.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""

from __future__ import print_function
from sys import stdout, maxsize
import csv

maxInt = 2147483647
decrement = True

while decrement:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.
    # http://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072

    decrement = False
    try:
        csv.field_size_limit(maxInt)
    except OverflowError:
        maxInt = int(maxInt/10)
        decrement = True

verb = False

"""
Doc :
    A class to read the SWATH scoring logs from Peakview and OpenSWATH

    The data is returned from the "parse_files" functions as a list of runs.
    Each run contains a number of precursors and each precursors contains
    peakgroups (e.g. peakgroups that were found in chromatographic space).

    usage: 
    reader = SWATHScoringReader.newReader(options.infiles, options.file_format)
    this_exp.runs = reader.parse_files(options.realign_runs)
"""

from msproteomicstoolslib.data_structures.PeakGroup import MinimalPeakGroup, GuiPeakGroup, GeneralPeakGroup
from msproteomicstoolslib.data_structures.Precursor import GeneralPrecursor, Precursor
from msproteomicstoolslib.data_structures.Run import Run
from msproteomicstoolslib.format.SWATHScoringMapper import inferMapping

class ReadFilter(object):
    """
    A callable class which can pre-filters a row and determine whether the row can be skipped.

    If the call returns true, the row is examined but if it returns false, the row should be skipped.
    """

    def __call__(self, row, header):
        return True
    def getSQL(self):
        return ""

#
# The Readers of the Scoring files
#

class SWATHScoringReader:

    def __init__(self):
        raise Exception("Abstract class")

    def parse_row(self, run, this_row, read_exp_RT):
        raise Exception("Abstract method")

    @staticmethod
    def newReader(infiles, filetype, readmethod="minimal",
                  readfilter=ReadFilter(), errorHandling="strict", enable_isotopic_grouping=False, read_cluster_id=True):
        """
        newReader(infiles, filetype, readmethod="minimal", readfilter=ReadFilter(), errorHandling="strict", enable_isotopic_grouping=False)

        Factory to create a new reader
        
        """
        if filetype  == "openswath": 
            if len(infiles) > 0 and infiles[0].endswith(".osw"):
                return OpenSWATH_OSW_SWATHScoringReader(infiles, readmethod,
                                                    readfilter, errorHandling,
                                                    enable_isotopic_grouping=enable_isotopic_grouping,
                                                    read_cluster_id=read_cluster_id)
            else:
                return OpenSWATH_SWATHScoringReader(infiles, readmethod,
                                                    readfilter, errorHandling,
                                                    enable_isotopic_grouping=enable_isotopic_grouping,
                                                    read_cluster_id=read_cluster_id)
        elif filetype  == "mprophet": 
            return mProphet_SWATHScoringReader(infiles, readmethod, readfilter)
        elif filetype  == "peakview": 
            return Peakview_SWATHScoringReader(infiles, readmethod, readfilter)
        elif filetype  == "peakview_preprocess": 
            return PeakviewPP_SWATHScoringReader(infiles, readmethod, readfilter)
        else:
            raise Exception("Unknown filetype '%s', allowed types are %s" % (decoy, str(filetypes) ) )

    def parse_files(self, read_exp_RT=True, verbosity=10, useCython=False):
        """Parse the input file(s) csv and osw.

        Args:
            read_exp_RT(bool) : to read the real, experimental retention time
                (default behavior) or the delta iRT should be used instead.
        
        Returns:
            runs(list(SWATHScoringReader.Run))

        A single CSV file might contain more than one run and thus to create
        unique run ids, we number the runs as xx_yy where xx is the current file
        number and yy is the run found in the current file. However, if an
        alignment has already been performed and each run has already obtained a
        unique run id, we can directly use the previous alignment id.
        """

        print("Parsing input files")
        from sys import stdout
        import csv
        skipped = 0; read = 0
        runs = []
        for file_nr, f in enumerate(self.infiles):
            if verbosity >= 10:
                stdout.write("\rReading %s" % str(f))
                stdout.flush()
            if f.endswith(".osw"):
                if not self.infiles_chromfiles_map:
                    read += self.parse_file(f, runs, useCython)
                else:
                    read += self.parse_file_selective(f, runs, useCython)
                continue

            header_dict = {}
            if f.endswith('.gz'):
                import gzip 
                filehandler = gzip.open(f,'rb')
            else:
                filehandler = open(f)
            reader = csv.reader(filehandler, delimiter="\t")
            header = next(reader)
            for i,n in enumerate(header):
                header_dict[n] = i
            if verbosity >= 10:
                stdout.write("\rReading file %s" % (str(f)) )
                stdout.flush()

            # Check if runs are already aligned (only one input file and correct header)
            already_aligned = (len(self.infiles) == 1 and self.aligned_run_id_name in header_dict)

            for this_row in reader:
                if already_aligned:
                    runid = this_row[header_dict[self.aligned_run_id_name]]
                else:
                    runnr = this_row[header_dict[self.run_id_name]]
                    runid = runnr + "_" + str(file_nr)

                current_run = [r for r in runs if r.get_id() == runid]
                # check if we have a new run
                if len(current_run) == 0:
                    orig_fname = None
                    aligned_fname = None
                    if "align_origfilename" in header_dict:
                        aligned_fname = this_row[header_dict[ "align_origfilename"] ]
                    if "filename" in header_dict:
                        orig_fname = this_row[header_dict[ "filename"] ]
                    current_run = Run(header, header_dict, runid, f, orig_fname, aligned_fname, useCython=useCython)
                    runs.append(current_run)
                    print(current_run, "maps to", orig_fname)
                else: 
                    assert len(current_run) == 1
                    current_run = current_run[0]

                if not self.readfilter(this_row, current_run.header_dict):
                    skipped += 1
                    continue

                read += 1
                # Unfortunately, since we are using csv, tell() will not work...
                # print("parse row at", filehandler.tell())
                self.parse_row(current_run, this_row, read_exp_RT)

        # Here we check that each run indeed has a unique id
        assert len(set([r.get_id() for r in runs])) == len(runs) # each run has a unique id
        if verbosity >= 10: stdout.write("\r\r\n") # clean up
        print("Found %s runs, read %s lines and skipped %s lines" % (len(runs), read, skipped))
        return runs
    
    def map_infiles_chromfiles(self, chromatogramFiles, useCython=False):
        """
        Updates the infiles_chromfiles_map with a dictionary. The dictionary has a mapping between each feature
        file (osw) and chromatogram file (chrom.mzML).
        """
        self.infiles_chromfiles_map = getMapping(chromatogramFiles, self.infiles, useCython)

class OpenSWATH_OSW_SWATHScoringReader(SWATHScoringReader):
    """
    Parser for OpenSWATH output in OSW (sqlite) format
    """

    def __init__(self, infiles, readmethod="minimal", readfilter=ReadFilter(), errorHandling="strict", enable_isotopic_grouping=False, read_cluster_id=True):
        self.infiles = infiles
        self.run_id_name = "run_id"
        self.readmethod = readmethod
        self.aligned_run_id_name = "align_runid"
        self.readfilter = readfilter
        self.errorHandling = errorHandling
        self.sequence_col = "Sequence"
        self.read_cluster_id = read_cluster_id
        self.infiles_chromfiles_map = {}
        if readmethod == "cminimal":
            try:
                from msproteomicstoolslib.cython.Precursor import CyPrecursor
                from msproteomicstoolslib.cython._optimized import CyPrecursorWrapperOnly
                self.Precursor = CyPrecursor
                self.Precursor = CyPrecursorWrapperOnly
            except ImportError as e:
                print ("Requested method 'cminimal' but Cython extensions seem to be missing. Please compile and add them or use readmethod 'minimal'")
                raise ValueError("Need Cython extensions for 'cminimal' readmethod.")

        elif readmethod == "minimal":
            self.Precursor = Precursor
        elif readmethod == "gui":
            self.Precursor = GeneralPrecursor
            self.PeakGroup = GuiPeakGroup
            self.sequence_col = "FullPeptideName"
        else:
            # complete
            self.Precursor = GeneralPrecursor
            self.PeakGroup = GeneralPeakGroup

        if enable_isotopic_grouping:
            self.peptide_group_label_name = "peptide_group_label"
        else:
            self.peptide_group_label_name = "transition_group_id"

    def parse_file(self, filename, runs, useCython):
        """
        Parse a whole OSW file (which may contain data from multiple runs)
        """

        import sqlite3
        conn = sqlite3.connect(filename)
        c = conn.cursor()

        # Retrieve and then iterate over all available runs
        query = """SELECT ID, FILENAME FROM RUN"""
        res = [row for row in c.execute(query)]
        nrows = 0
        for row in res:
            runid = row[0]
            current_run = Run([], {}, runid, filename, row[1], row[1], useCython=useCython)
            runs.append(current_run)
            nrows += self._parse_file(filename, current_run, runid, conn)
        return nrows

    def parse_file_selective(self, filename, runs, useCython):
        """
        Parse a OSW file. Pick runs that has corresponding mzML file. 
        """

        nrows = 0
        if filename not in self.infiles_chromfiles_map:
            return nrows

        import sqlite3
        conn = sqlite3.connect(filename)

        # Retrieve available runs from the dictionary and then iterate over them
        for current_run in self.infiles_chromfiles_map[filename]:
            runid = current_run.get_id()
            runs.append(current_run)
            nrows += self._parse_file(filename, current_run, runid, conn)
        conn.close()
        return nrows

    def _parse_file(self, filename, run, run_id, conn):

        # Make sure the correct indices exist
        c = conn.cursor()
        c.executescript('''
            CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
            CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
            CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
            CREATE INDEX IF NOT EXISTS idx_feature_run_id ON FEATURE (RUN_ID);
            CREATE INDEX IF NOT EXISTS idx_feature_ms2_feature_id ON FEATURE_MS2 (FEATURE_ID);
            CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);
        ''')
        # TODO: Shubham: Add left width and right width
        query = """
        SELECT PRECURSOR.ID, SCORE_MS2.QVALUE, FEATURE.EXP_RT, PRECURSOR.DECOY, PEPTIDE.MODIFIED_SEQUENCE, PEPTIDE.ID, FEATURE.ID, AREA_INTENSITY, SCORE_MS2.SCORE
        FROM SCORE_MS2
        INNER JOIN (SELECT FEATURE_ID, AREA_INTENSITY FROM FEATURE_MS2) AS FEATURE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE_MS2.FEATURE_ID
        INNER JOIN (SELECT ID, PRECURSOR_ID, EXP_RT, RUN_ID FROM FEATURE) AS FEATURE ON SCORE_MS2.FEATURE_ID = FEATURE.ID
        INNER JOIN (SELECT ID, DECOY FROM PRECURSOR) AS PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
        INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID
        INNER JOIN (SELECT ID, MODIFIED_SEQUENCE FROM PEPTIDE) AS PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
        WHERE RUN_ID = %s
        %s
        """ % (run_id, self.readfilter.getSQL() )
        q = [row for row in c.execute(query)]
        for row in q:

            trgr_id = row[0]
            unique_peakgroup_id = row[0]
            sequence = row[4]
            peptide_group_label = row[5]

            # Attributes that only need to be present in strict mode
            intensity = row[7]
            d_score = row[8]
            cluster_id = -1
            retention_time = row[2]
            fdr_score = row[1]
            protein_name = "NA"
            thisid = row[6]
            if row[3] == 0:
                decoy = "FALSE"
            elif row[3] == 1:
                decoy = "TRUE"
            else:
                raise Exception("Unknown decoy", row[3])
            # left_width = run [8]
            # right_width = run [9]

            # If the peptide does not yet exist, generate it
            if not run.hasPrecursor(peptide_group_label, trgr_id):
                p = self.Precursor(trgr_id, run)
                p.setProteinName(protein_name)
                p.setSequence(sequence)
                p.set_decoy(decoy)
                run.addPrecursor(p, peptide_group_label)

            if self.readmethod == "cminimal":
                peakgroup_tuple = (thisid, fdr_score, retention_time, intensity, d_score)
                run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id, cluster_id)
            elif self.readmethod == "minimal":
                peakgroup_tuple = (thisid, fdr_score, retention_time, intensity, d_score)
                run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id, cluster_id)
            elif self.readmethod == "gui":
                leftWidth = this_row[run.header_dict[left_width_name]]
                assay_rt = -1
                if "assay_rt" in run.header_dict:
                    assay_rt = this_row[run.header_dict["assay_rt"]]
                rightWidth = this_row[run.header_dict[right_width_name]]
                charge = this_row[run.header_dict[charge_name]]
                peakgroup = self.PeakGroup(fdr_score, intensity, leftWidth, rightWidth, assay_rt, run.getPrecursor(peptide_group_label, trgr_id))
                peakgroup.charge = charge
                run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)
            elif self.readmethod == "complete":
                peakgroup = self.PeakGroup(this_row, run, run.getPrecursor(peptide_group_label, trgr_id))
                peakgroup.set_normalized_retentiontime(retention_time)
                peakgroup.set_fdr_score(fdr_score)
                peakgroup.set_feature_id(thisid)
                peakgroup.set_intensity(intensity)
                peakgroup.setClusterID(cluster_id)
                run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)
            else:
                raise Exception("Unknown readmethod", self.readmethod)
        return len(q)
            
class OpenSWATH_SWATHScoringReader(SWATHScoringReader):
    """
    Parser for OpenSWATH output
    """

    def __init__(self, infiles, readmethod="minimal", readfilter=ReadFilter(), errorHandling="strict", enable_isotopic_grouping=False, read_cluster_id=True):
        self.infiles = infiles
        self.run_id_name = "run_id"
        self.readmethod = readmethod
        self.aligned_run_id_name = "align_runid"
        self.readfilter = readfilter
        self.errorHandling = errorHandling
        self.sequence_col = "Sequence"
        self.read_cluster_id = read_cluster_id
        if readmethod == "cminimal":
            try:
                from msproteomicstoolslib.cython.data_structures import CyPrecursorWrapperOnly
                self.Precursor = CyPrecursorWrapperOnly
            except ImportError as e:
                print ("Requested method 'cminimal' but Cython extensions seem to be missing. Please compile and add them or use readmethod 'minimal'")
                raise ValueError("Need Cython extensions for 'cminimal' readmethod.")

        elif readmethod == "minimal":
            self.Precursor = Precursor
        elif readmethod == "gui":
            self.Precursor = GeneralPrecursor
            self.PeakGroup = GuiPeakGroup
            self.sequence_col = "FullPeptideName"
        else:
            # complete
            self.Precursor = GeneralPrecursor
            self.PeakGroup = GeneralPeakGroup

        if enable_isotopic_grouping:
            self.peptide_group_label_name = "peptide_group_label"
        else:
            self.peptide_group_label_name = "transition_group_id"

    def parse_row(self, run, this_row, read_exp_RT):
        decoy_name = "decoy"
        fdr_score_name = "m_score"
        fdr_score_name_alt = "transition_group_id_m_score"
        fdr_score_name_alt2 = "q_value"
        dscore_name = "d_score"
        unique_peakgroup_id_name = "transition_group_id"
        rt_name = "delta_rt"
        run_id_name = "run_id"
        protein_id_col = "ProteinName"
        unique_feature_id_name = "id"
        intensity_name = "Intensity"
        decoy = "FALSE"
        left_width_name = "leftWidth"
        right_width_name = "rightWidth"
        charge_name = "Charge"
        cluster_id = -1

        # use the aligned retention time if it is available!
        if "aligned_rt" in run.header_dict: 
            rt_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            rt_name = "RT"
        if "align_clusterid" in run.header_dict and self.read_cluster_id:
            cluster_id = int(this_row[run.header_dict["align_clusterid"]])

        trgr_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        unique_peakgroup_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        sequence = this_row[run.header_dict[self.sequence_col]]
        peptide_group_label = trgr_id

        if self.peptide_group_label_name in run.header_dict: 
            peptide_group_label = this_row[run.header_dict[self.peptide_group_label_name]]

        # Attributes that only need to be present in strict mode
        retention_time = -1
        fdr_score = -1
        protein_name = "NA"
        thisid = -1
        
        # stop errors occuring if the value passed is empty
        to_float = lambda x: float('nan' if x == '' else x)
        
        if fdr_score_name in run.header_dict:
            fdr_score = to_float(this_row[run.header_dict[fdr_score_name]])
        elif fdr_score_name_alt in run.header_dict:
            fdr_score = to_float(this_row[run.header_dict[fdr_score_name_alt]])
        elif fdr_score_name_alt2 in run.header_dict:
            fdr_score = to_float(this_row[run.header_dict[fdr_score_name_alt2]])
        elif self.errorHandling == "strict": 
            raise Exception("Did not find essential column " + fdr_score_name + " or " + fdr_score_name_alt + " or " + fdr_score_name_alt2)

        try:
            thisid = this_row[run.header_dict[unique_feature_id_name]]
            retention_time = to_float(this_row[run.header_dict[rt_name]])
            d_score = to_float(this_row[run.header_dict[dscore_name]])
            protein_name = this_row[run.header_dict[protein_id_col]]
        except KeyError:
            if self.errorHandling == "strict": 
                raise Exception("Did not find essential column.")

        # Optional attributes
        intensity = -1
        if intensity_name in run.header_dict:
            intensity = to_float(this_row[run.header_dict[intensity_name]])
        if "decoy" in run.header_dict:
            decoy = this_row[run.header_dict[decoy_name]]

        # If the peptide does not yet exist, generate it
        if not run.hasPrecursor(peptide_group_label, trgr_id):
            p = self.Precursor(trgr_id, run)
            p.setProteinName(protein_name)
            p.setSequence(sequence)
            p.set_decoy(decoy)
            run.addPrecursor(p, peptide_group_label)

        if self.readmethod == "cminimal":
            peakgroup_tuple = (thisid, fdr_score, retention_time, intensity, d_score)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id, cluster_id)
        elif self.readmethod == "minimal":
            peakgroup_tuple = (thisid, fdr_score, retention_time, intensity, d_score)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id, cluster_id)
        elif self.readmethod == "gui":
            leftWidth = this_row[run.header_dict[left_width_name]]
            assay_rt = -1
            if "assay_rt" in run.header_dict:
                assay_rt = this_row[run.header_dict["assay_rt"]]
            rightWidth = this_row[run.header_dict[right_width_name]]
            charge = this_row[run.header_dict[charge_name]]
            peakgroup = self.PeakGroup(fdr_score, intensity, leftWidth, rightWidth, assay_rt, run.getPrecursor(peptide_group_label, trgr_id))
            peakgroup.charge = charge
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)
        elif self.readmethod == "complete":
            peakgroup = self.PeakGroup(this_row, run, run.getPrecursor(peptide_group_label, trgr_id))
            peakgroup.set_normalized_retentiontime(retention_time)
            peakgroup.set_fdr_score(fdr_score)
            peakgroup.set_feature_id(thisid)
            peakgroup.set_intensity(intensity)
            peakgroup.setClusterID(cluster_id)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)
        else:
            raise Exception("Unknown readmethod", self.readmethod)

class mProphet_SWATHScoringReader(SWATHScoringReader):
    """
    Parser for mProphet output
    """

    def __init__(self, infiles, readmethod="minimal", readfilter=ReadFilter(), enable_isotopic_grouping=False):
        self.infiles = infiles
        self.run_id_name = "run_id"
        self.readmethod = readmethod
        self.aligned_run_id_name = "align_runid"
        self.readfilter = readfilter
        if readmethod == "minimal":
            self.Precursor = Precursor
        else:
            self.Precursor = GeneralPrecursor
            self.PeakGroup = GeneralPeakGroup

        if enable_isotopic_grouping:
            raise Exception("Cannot use isotopic grouping with mProphet data.")

    def parse_row(self, run, this_row, read_exp_RT):
        decoy_name = "decoy"
        fdr_score_name = "m_score"
        unique_peakgroup_id_name = "transition_group_id"
        # rt_name = "delta_rt"
        run_id_name = "run_id"
        protein_id_col = "protein"
        sequence_col = "transition_group_pepseq"
        intensity_name = "log10_max_apex_intensity"
        decoy = "FALSE"

        # use the aligned retention time if it is available!
        if "aligned_rt" in run.header_dict: 
            rt_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            rt_name = "Tr" 
            retention_time = float(this_row[run.header_dict["Tr"]]) 
        else:
            retention_time = float(this_row[run.header_dict["iRT_empirical"]]) - float(this_row[run.header_dict["iRT_predicted"]])

        # create some id
        import uuid 
        thisid = str(uuid.uuid1() )
        # thisid = this_row[run.header_dict[unique_feature_id_name]]

        trgr_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        protein_name = this_row[run.header_dict[protein_id_col]]
        sequence = this_row[run.header_dict[sequence_col]]
        fdr_score = float(this_row[run.header_dict[fdr_score_name]])
        unique_peakgroup_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        intensity = -1
        if intensity_name in run.header_dict:
            intensity = float(this_row[run.header_dict[intensity_name]])
        if "decoy" in run.header_dict:
            decoy = this_row[run.header_dict[decoy_name]]
        run_id = this_row[run.header_dict[run_id_name]]

        # If the peptide does not yet exist, generate it
        peptide_group_label = trgr_id
        if not run.hasPrecursor(peptide_group_label, trgr_id):
            p = self.Precursor(trgr_id, run)
            p.protein_name = protein_name
            p.sequence = sequence
            p.set_decoy(decoy)
            run.addPrecursor(p, peptide_group_label)
        if self.readmethod == "minimal":
            peakgroup_tuple = (thisid, fdr_score, retention_time, intensity)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id)
        else:
            peakgroup = self.PeakGroup(this_row, run, run.getPrecursor(peptide_group_label, trgr_id))
            peakgroup.set_normalized_retentiontime(retention_time)
            peakgroup.set_fdr_score(fdr_score)
            peakgroup.set_feature_id(thisid)
            peakgroup.set_intensity(intensity)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)

class Peakview_SWATHScoringReader(SWATHScoringReader):
    """
    Parser for Peakview output
    """

    def __init__(self, infiles, readmethod="minimal", readfilter=ReadFilter(), enable_isotopic_grouping=False):
        self.infiles = infiles
        self.run_id_name = "Sample"
        self.aligned_run_id_name = "align_runid"
        self.readmethod = readmethod
        self.readfilter = readfilter

        # Only minimal reading is implemented
        if readmethod == "minimal":
            self.Precursor = Precursor
        else:
            raise Exception("Cannot only use readmethod minimal Peakview data.")

        if enable_isotopic_grouping:
            raise Exception("Cannot use isotopic grouping with Peakview data.")


    def parse_row(self, run, this_row, read_exp_RT):
        decoy_name = "Decoy"
        fdr_score_name = "Score" # TODO invert score!!!
        unique_peakgroup_id_name = "Peptide" ## Does not exist!!
        run_id_name = "Sample"
        protein_id_col = "Protein"
        sequence_col = "Peptide"
        unique_feature_id_name = "id" # does not exist!!!
        decoy = "FALSE"
        intensity_name = "MaxPeak.Intensity"

        rt_name = "empirical_iRT"
        if not rt_name in run.header_dict:
            rt_name = "Median RT"

        # use the aligned retention time if it is available!
        if "aligned_rt" in run.header_dict: 
            rt_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            rt_name = "Median RT"

        # create some id
        import uuid 
        thisid = str(uuid.uuid1() )
        # thisid = this_row[run.header_dict[unique_feature_id_name]]

        if len(this_row) < run.header_dict[fdr_score_name] :
            # what to do here!? 
            return

        protein_name = this_row[run.header_dict[protein_id_col]]
        sequence = this_row[run.header_dict[sequence_col]]
        trgr_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        unique_peakgroup_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        intensity = float(this_row[run.header_dict[intensity_name]])
        #print run.header_dict
        #print run.header_dict["Score"]

        # compute 1/score to have a score that gets better when smaller
        fdr_score = float(this_row[run.header_dict[fdr_score_name]])
        fdr_score = 1/fdr_score

        retention_time = float(this_row[run.header_dict[rt_name]])
        if "decoy" in run.header_dict:
            decoy = this_row[run.header_dict[decoy_name]]
        run_id = this_row[run.header_dict[run_id_name]]

        # If the peptide does not yet exist, generate it
        peptide_group_label = trgr_id
        if not run.hasPrecursor(peptide_group_label, trgr_id):
            p = self.Precursor(trgr_id, run)
            p.protein_name = protein_name
            p.sequence = sequence
            p.set_decoy(decoy)
            run.addPrecursor(p, peptide_group_label)
            if verb: print("add peptide", trgr_id)

        # Only minimal reading is implemented
        if self.readmethod == "minimal":
            if verb: print("append tuple", peakgroup_tuple)
            peakgroup_tuple = (thisid, fdr_score, retention_time,intensity)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id)
        else:
            raise NotImplemented

class PeakviewPP_SWATHScoringReader(Peakview_SWATHScoringReader):
    """
    Parser for Peakview output
    """

    def __init__(self, infiles, readmethod="minimal", readfilter=ReadFilter(), enable_isotopic_grouping=False):
        self.infiles = infiles
        self.run_id_name = "Sample"
        self.aligned_run_id_name = "align_runid"
        self.readmethod = readmethod
        self.readfilter = readfilter
        self.errorHandling = "strict"

        # Only minimal reading is implemented
        if readmethod == "minimal":
            self.Precursor = Precursor
        else:
            raise Exception("Cannot only use readmethod minimal Peakview data.")

        if enable_isotopic_grouping:
            raise Exception("Cannot use isotopic grouping with Peakview data.")
        else:
            self.peptide_group_label_name = "Pep Index"

    def parse_row(self, run, this_row, read_exp_RT):
        decoy_name = "decoy"
        fdr_score_name = "m_score"
        ## dscore_name = "d_score"
        unique_peakgroup_id_name = "transition_group_id"
        ## rt_name = "delta_rt"
        run_id_name = "run_id"
        protein_id_col = "ProteinName"
        intensity_name = "Intensity"
        decoy = "FALSE"
        # left_width_name = "leftWidth"
        # right_width_name = "rightWidth"
        charge_name = "Charge"
        cluster_id = -1

        decoy_name = "Decoy"
        fdr_score_name = "Score" # note this score is better when higher
        unique_peakgroup_id_name = "Pep Index"
        rt_name = "Median RT"
        run_id_name = "Sample"
        protein_id_col = "Protein"
        self.sequence_col = "Peptide"
        unique_feature_id_name = "preprocess_id"
        decoy = "FALSE"
        intensity_name = "MaxPeak Intensity"
        charge_name = "Precursor Charge"
        d_score = 2

        # use the aligned retention time if it is available!
        if "aligned_rt" in run.header_dict: 
            rt_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            rt_name = "Median RT"
        if "align_clusterid" in run.header_dict: 
            cluster_id = int(this_row[run.header_dict["align_clusterid"]])

        trgr_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        unique_peakgroup_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        sequence = this_row[run.header_dict[self.sequence_col]]
        peptide_group_label = trgr_id

        if self.peptide_group_label_name in run.header_dict: 
            peptide_group_label = this_row[run.header_dict[self.peptide_group_label_name]]

        # Attributes that only need to be present in strict mode
        retention_time = -1
        fdr_score = -1
        protein_name = "NA"
        thisid = -1
        try:
            fdr_score = float(this_row[run.header_dict[fdr_score_name]])
            protein_name = this_row[run.header_dict[protein_id_col]]
            thisid = this_row[run.header_dict[unique_feature_id_name]]
            retention_time = float(this_row[run.header_dict[rt_name]])
        except KeyError:
            if self.errorHandling == "strict": 
                raise Exception("Did not find essential column.")

        # compute 1/score to have a score that gets better when smaller
        fdr_score = 1/fdr_score

        # Optional attributes
        intensity = -1
        if intensity_name in run.header_dict:
            intensity = float(this_row[run.header_dict[intensity_name]])
        if "Decoy" in run.header_dict:
            decoy = this_row[run.header_dict[decoy_name]]
            if decoy in ["True", "TRUE"]:
                decoy = "TRUE"
            elif decoy in ["False", "FALSE"]:
                decoy = "FALSE"

        # If the peptide does not yet exist, generate it
        if not run.hasPrecursor(peptide_group_label, trgr_id):
            p = self.Precursor(trgr_id, run)
            p.protein_name = protein_name
            p.sequence = sequence
            p.set_decoy(decoy)
            run.addPrecursor(p, peptide_group_label)

        if self.readmethod == "minimal":
            peakgroup_tuple = (thisid, fdr_score, retention_time, intensity, d_score)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id, cluster_id)
        elif self.readmethod == "gui":
            leftWidth = this_row[run.header_dict[left_width_name]]
            rightWidth = this_row[run.header_dict[right_width_name]]
            charge = this_row[run.header_dict[charge_name]]
            peakgroup = self.PeakGroup(fdr_score, intensity, leftWidth, rightWidth, run.getPrecursor(peptide_group_label, trgr_id))
            peakgroup.charge = charge
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)
        elif self.readmethod == "complete":
            peakgroup = self.PeakGroup(this_row, run, run.getPrecursor(peptide_group_label, trgr_id))
            peakgroup.set_normalized_retentiontime(retention_time)
            peakgroup.set_fdr_score(fdr_score)
            peakgroup.set_feature_id(thisid)
            peakgroup.set_intensity(intensity)
            peakgroup.setClusterID(cluster_id)
            run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)

def getMapping(chromatogramFiles, featureFiles, useCython = False):
    """
    Returns a dictionary of feature files with Run objects. Only those Runs are included which
    have corresponding mzML files available.
    Sometimes, a single feature file is provided for all chromatogram files, this function checks
    whether there is corresponding chromatogram file present for each run.
    
    >>> chromatogramFiles = ["file1.chrom.mzML", "file2.chrom.mzML"]
    >>> featureFiles = ["file1.osw", "file2.osw"]
    >>> featureFiles_chromFiles_map = getMapping(chromatogramFiles, featureFiles)
    >>> {"file1.osw" : [run1], "file2.osw" : [run2]}
    """

    # Read osw file. Create a dictionary of "Feature file": Run()
    MS_feature_fileMapping = getRunfromFeatureFile(featureFiles, useCython)

    # Read mzML, get a list of all provided mzML file.
    chromFiles = []
    for filename in chromatogramFiles:
        base_name = getBaseName(filename)
        if base_name in chromFiles:
            print(str(filename) + " has basename same as of another mzML file. Basename is obtained by removing directory separator, .gz, .sqMass, .mzML and .chrom extensions.")
        else:
            chromFiles.append(base_name)

    # Iterate through the osw file and find associated mzML file.
    # Remove Run from the dictionary if mzML file is missing and Throw warning.
    # Add to the final dictionary if the basename has associated chromatogram file and OpenSWATH feature file.

    for featureFile, runs in MS_feature_fileMapping.items():
        for run in runs:
            MSrun = getBaseName(run.get_openswath_filename())
            if(MSrun not in chromFiles):
                # Warning
                print("Chromatogram file associated with " + run.get_openswath_filename() + " from " + run.get_original_filename() + " is not found. Skipping!")
                MS_feature_fileMapping[featureFile].remove(run)
        if len(MS_feature_fileMapping[featureFile]) == 0:
            del MS_feature_fileMapping[featureFile]
    
    if not MS_feature_fileMapping:
        # Error
        raise Exception("Feature files and chromatogram files do not have matching names.")
    return MS_feature_fileMapping

def getRunfromFeatureFile(featureFiles, useCython = False):
    """
    Return as dictionary with key as feature file and value as associated Run objects.

    >>> featureFiles = ["merged.osw"]
    >>> fileMapping = getRunfromFeatureFile(featureFiles)
    >>> fileMapping = {"merged.osw": [Run0, Run1, Run2]}
    """

    import sqlite3
    MSfile_featureFile_mapping = {}
    for filename in featureFiles:
        MSfile_featureFile_mapping[filename] = []
        conn = sqlite3.connect(filename)
        c = conn.cursor()
        try:
            # Retrieve and then iterate over all available runs
            query = """SELECT ID, FILENAME FROM RUN"""
            results = [row for row in c.execute(query)]
            for (run_id, MS_file) in results:
                current_run = Run([], {}, run_id, filename, MS_file, MS_file, useCython=useCython)
                MSfile_featureFile_mapping[filename].append(current_run)
        
        except sqlite3.Error as e:
            print("An error occured in reading file " + str(filename) + ", ", e.args[0])
            conn.close()
                
        # Close the connection
        conn.close()
    return MSfile_featureFile_mapping

def getBaseName(filename):
    """
    Returns a basename for filename. Basename is obtained by removing directory separator, .gz, .sqMass, .mzML and .chrom extensions.
    
    >>> getBaseName('data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz')
    >>> hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt
    """

    fileBase = filename
    # Remove directory separator from filename
    for dir_sep in ["/", "\\"]:
        fileBase = fileBase.split(dir_sep)[-1]
    # First remove gz, then sqmass, then mzML & then chrom, if present.
    for ending in [".gz", ".sqMass", ".mzML", ".chrom"]:
        fileBase = fileBase.split(ending)[0]
    return fileBase

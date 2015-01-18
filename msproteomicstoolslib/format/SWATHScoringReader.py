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

from sys import stdout
import csv

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

class ReadFilter(object):
    """
    A callable class which can pre-filters a row and determine whether the row can be skipped.

    If the call returns true, the row is examined but if it returns false, the row should be skipped.
    """

    def __call__(self, row, header):
        return True

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
                  readfilter=ReadFilter(), errorHandling="strict", enable_isotopic_grouping=False):
        """
        newReader(infiles, filetype, readmethod="minimal", readfilter=ReadFilter(), errorHandling="strict", enable_isotopic_grouping=False)

        Factory to create a new reader
        
        """
        if filetype  == "openswath": 
            return OpenSWATH_SWATHScoringReader(infiles, readmethod,
                                                readfilter, errorHandling,
                                                enable_isotopic_grouping=enable_isotopic_grouping)
        elif filetype  == "mprophet": 
            return mProphet_SWATHScoringReader(infiles, readmethod, readfilter)
        elif filetype  == "peakview": 
            return Peakview_SWATHScoringReader(infiles, readmethod, readfilter)
        elif filetype  == "peakview_preprocess": 
            return PeakviewPP_SWATHScoringReader(infiles, readmethod, readfilter)
        else:
            raise Exception("Unknown filetype '%s', allowed types are %s" % (decoy, str(filetypes) ) )

    def parse_files(self, read_exp_RT=True, verbosity=10):
      """Parse the input file(s) (CSV).

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

      print "Parsing input files"
      from sys import stdout
      import csv
      skipped = 0; read = 0
      runs = []
      for file_nr, f in enumerate(self.infiles):
        if verbosity >= 10:
            stdout.write("\rReading %s" % str(f))
            stdout.flush()
        header_dict = {}
        if f.endswith('.gz'):
            import gzip 
            filehandler = gzip.open(f,'rb')
        else:
            filehandler = open(f)
        reader = csv.reader(filehandler, delimiter="\t")
        header = reader.next()
        for i,n in enumerate(header):
          header_dict[n] = i
        if verbosity >= 10:
            stdout.write("\rReading file %s" % (str(f)) )
            stdout.flush()

        # Check if runs are already aligned (only one input file and correct header)
        already_aligned = (len(self.infiles) == 1 and header_dict.has_key(self.aligned_run_id_name))

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
                if header_dict.has_key("align_origfilename"):
                    aligned_fname = this_row[header_dict[ "align_origfilename"] ]
                if header_dict.has_key("filename"):
                    orig_fname = this_row[header_dict[ "filename"] ]
                current_run = Run(header, header_dict, runid, f, orig_fname, aligned_fname)
                runs.append(current_run)
                print current_run, "maps to ", orig_fname
            else: 
                assert len(current_run) == 1
                current_run = current_run[0]

            if not self.readfilter(this_row, current_run.header_dict):
                skipped += 1
                continue

            read += 1
            # Unfortunately, since we are using csv, tell() will not work...
            # print "parse row at", filehandler.tell()
            self.parse_row(current_run, this_row, read_exp_RT)

      # Here we check that each run indeed has a unique id
      assert len(set([r.get_id() for r in runs])) == len(runs) # each run has a unique id
      if verbosity >= 10: stdout.write("\r\r\n") # clean up
      print "Found %s runs, read %s lines and skipped %s lines" % (len(runs), read, skipped)
      return runs

class OpenSWATH_SWATHScoringReader(SWATHScoringReader):
    """
    Parser for OpenSWATH output
    """

    def __init__(self, infiles, readmethod="minimal", readfilter=ReadFilter(), errorHandling="strict", enable_isotopic_grouping=False):
        self.infiles = infiles
        self.run_id_name = "run_id"
        self.readmethod = readmethod
        self.aligned_run_id_name = "align_runid"
        self.readfilter = readfilter
        self.errorHandling = errorHandling
        self.sequence_col = "Sequence"
        if readmethod == "minimal":
            self.Precursor = Precursor
        elif readmethod == "gui":
            self.Precursor = GeneralPrecursor
            self.PeakGroup = GuiPeakGroup
            self.sequence_col = "FullPeptideName"
        else:
            self.Precursor = GeneralPrecursor
            self.PeakGroup = GeneralPeakGroup

        if enable_isotopic_grouping:
            self.peptide_group_label_name = "peptide_group_label"
        else:
            self.peptide_group_label_name = "transition_group_id"

    def parse_row(self, run, this_row, read_exp_RT):
        decoy_name = "decoy"
        fdr_score_name = "m_score"
        dscore_name = "d_score"
        unique_peakgroup_id_name = "transition_group_id"
        diff_from_assay_in_sec_name = "delta_rt"
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
            diff_from_assay_in_sec_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            diff_from_assay_in_sec_name = "RT"
        if "align_clusterid" in run.header_dict: 
            cluster_id = int(this_row[run.header_dict["align_clusterid"]])

        trgr_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        unique_peakgroup_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        sequence = this_row[run.header_dict[self.sequence_col]]
        peptide_group_label = trgr_id

        if self.peptide_group_label_name in run.header_dict: 
            peptide_group_label = this_row[run.header_dict[self.peptide_group_label_name]]

        # Attributes that only need to be present in strict mode
        diff_from_assay_seconds = -1
        fdr_score = -1
        protein_name = "NA"
        thisid = -1
        try:
            fdr_score = float(this_row[run.header_dict[fdr_score_name]])
            #fdr_score = 0.0001
            protein_name = this_row[run.header_dict[protein_id_col]]
            thisid = this_row[run.header_dict[unique_feature_id_name]]
            diff_from_assay_seconds = float(this_row[run.header_dict[diff_from_assay_in_sec_name]])
            d_score = float(this_row[run.header_dict[dscore_name]])
            #d_score = 2
        except KeyError:
            if self.errorHandling == "strict": 
                raise Exception("Did not find essential column.")

        # Optional attributes
        intensity = -1
        if run.header_dict.has_key(intensity_name):
            intensity = float(this_row[run.header_dict[intensity_name]])
        if "decoy" in run.header_dict:
            decoy = this_row[run.header_dict[decoy_name]]

        # If the peptide does not yet exist, generate it
        if not run.hasPrecursor(peptide_group_label, trgr_id):
          p = self.Precursor(trgr_id, run)
          p.protein_name = protein_name
          p.sequence = sequence
          p.set_decoy(decoy)
          run.addPrecursor(p, peptide_group_label)

        if self.readmethod == "minimal":
          peakgroup_tuple = (thisid, fdr_score, diff_from_assay_seconds, intensity, d_score)
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
          peakgroup.set_normalized_retentiontime(diff_from_assay_seconds)
          peakgroup.set_fdr_score(fdr_score)
          peakgroup.set_feature_id(thisid)
          peakgroup.set_intensity(intensity)
          peakgroup.setClusterID(cluster_id)
          run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)

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
        # diff_from_assay_in_sec_name = "delta_rt"
        run_id_name = "run_id"
        protein_id_col = "protein"
        sequence_col = "transition_group_pepseq"
        intensity_name = "log10_max_apex_intensity"
        decoy = "FALSE"

        # use the aligned retention time if it is available!
        if "aligned_rt" in run.header_dict: 
            diff_from_assay_in_sec_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            diff_from_assay_in_sec_name = "Tr" 
            diff_from_assay_seconds = float(this_row[run.header_dict["Tr"]]) 
        else:
            diff_from_assay_seconds = float(this_row[run.header_dict["iRT_empirical"]]) - float(this_row[run.header_dict["iRT_predicted"]])

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
        if run.header_dict.has_key(intensity_name):
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
          peakgroup_tuple = (thisid, fdr_score, diff_from_assay_seconds, intensity)
          run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id)
        else:
          peakgroup = self.PeakGroup(this_row, run, run.getPrecursor(peptide_group_label, trgr_id))
          peakgroup.set_normalized_retentiontime(diff_from_assay_seconds)
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

        diff_from_assay_in_sec_name = "empirical_iRT"
        if not diff_from_assay_in_sec_name in run.header_dict:
            diff_from_assay_in_sec_name = "Median RT"

        # use the aligned retention time if it is available!
        if "aligned_rt" in run.header_dict: 
            diff_from_assay_in_sec_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            diff_from_assay_in_sec_name = "Median RT"

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

        diff_from_assay_seconds = float(this_row[run.header_dict[diff_from_assay_in_sec_name]])
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
          if verb: print "add peptide", trgr_id

        # Only minimal reading is implemented
        if self.readmethod == "minimal":
          if verb: print "append tuple", peakgroup_tuple
          peakgroup_tuple = (thisid, fdr_score, diff_from_assay_seconds,intensity)
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
        ## diff_from_assay_in_sec_name = "delta_rt"
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
        diff_from_assay_in_sec_name = "Median RT"
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
            diff_from_assay_in_sec_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if read_exp_RT: 
            diff_from_assay_in_sec_name = "Median RT"
        if "align_clusterid" in run.header_dict: 
            cluster_id = int(this_row[run.header_dict["align_clusterid"]])

        trgr_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        unique_peakgroup_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        sequence = this_row[run.header_dict[self.sequence_col]]
        peptide_group_label = trgr_id

        if self.peptide_group_label_name in run.header_dict: 
            peptide_group_label = this_row[run.header_dict[self.peptide_group_label_name]]

        # Attributes that only need to be present in strict mode
        diff_from_assay_seconds = -1
        fdr_score = -1
        protein_name = "NA"
        thisid = -1
        try:
            fdr_score = float(this_row[run.header_dict[fdr_score_name]])
            protein_name = this_row[run.header_dict[protein_id_col]]
            thisid = this_row[run.header_dict[unique_feature_id_name]]
            diff_from_assay_seconds = float(this_row[run.header_dict[diff_from_assay_in_sec_name]])
        except KeyError:
            if self.errorHandling == "strict": 
                raise Exception("Did not find essential column.")

        # compute 1/score to have a score that gets better when smaller
        fdr_score = 1/fdr_score

        # Optional attributes
        intensity = -1
        if run.header_dict.has_key(intensity_name):
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
          peakgroup_tuple = (thisid, fdr_score, diff_from_assay_seconds, intensity, d_score)
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
          peakgroup.set_normalized_retentiontime(diff_from_assay_seconds)
          peakgroup.set_fdr_score(fdr_score)
          peakgroup.set_feature_id(thisid)
          peakgroup.set_intensity(intensity)
          peakgroup.setClusterID(cluster_id)
          run.getPrecursor(peptide_group_label, trgr_id).add_peakgroup(peakgroup)

def simpleInferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping,
                 sequences_mapping, protein_mapping, verbose=False):

    assert len(aligned_pg_files) == 1, "There should only be one file in simple mode"
    f = aligned_pg_files[0]

    # Produce simple mapping between runs and files (assume each file is one run)
    for i,raw in enumerate(rawdata_files):
        mapping[str(i)] = [ raw ]

    # Get the compression
    if f.endswith('.gz'):
        import gzip 
        filehandler = gzip.open(f,'rb')
    else:
        filehandler = open(f)

    # Get the dialect
    dialect = csv.Sniffer().sniff(filehandler.readline(), [',',';','\t'])
    filehandler.seek(0) 
    reader = csv.reader(filehandler, dialect)
    header = reader.next()

    header_dict = {}
    for i,n in enumerate(header):
        header_dict[n] = i

    for this_row in reader:
        peptide_name = this_row [ header_dict["chromatogram_super_group_id"]]
        precursor_name = this_row [ header_dict["chromatogram_group_id"]]
        transition_name = this_row [ header_dict["chromatogram_id"]]

        # Fill the sequence mapping
        tmp = sequences_mapping.get(peptide_name, [])
        if precursor_name not in tmp:
            tmp.append(precursor_name)
        sequences_mapping[peptide_name] = tmp

        # Fill the precursor mapping
        tmp = precursors_mapping.get(precursor_name, [])
        if transition_name not in tmp:
            tmp.append(transition_name)
        precursors_mapping[precursor_name] = tmp

def tramlInferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping,
                 sequences_mapping, protein_mapping, verbose=False):
    try:
        import pyopenms
    except ImportError as e:
        print "\nError!"
        print "Could not import pyOpenMS while trying to load a TraML file - please make sure pyOpenMS is installed."
        print "pyOpenMS is available from https://pypi.python.org/pypi/pyopenms"
        print 
        raise e

    assert len(aligned_pg_files) == 1, "There should only be one file in simple mode"
    f = aligned_pg_files[0]

    # Produce simple mapping between runs and files (assume each file is one run)
    for i,raw in enumerate(rawdata_files):
        mapping[str(i)] = [ raw ]

    targexp = pyopenms.TargetedExperiment()
    pyopenms.TraMLFile().load(f, targexp)

    for peptide_precursor in targexp.getPeptides():

        # Fill the protein mapping
        protein_id = peptide_precursor.protein_refs
        if len(protein_id) > 0:
            protein_id = protein_id[0] # take the first one ... 

        tmp = protein_mapping.get(protein_id, [])
        if peptide_precursor.sequence not in tmp:
            tmp.append(peptide_precursor.sequence)
        protein_mapping[protein_id] = tmp

        # Fill the sequence mapping
        tmp = sequences_mapping.get(peptide_precursor.sequence, [])
        if peptide_precursor.id not in tmp:
            tmp.append(peptide_precursor.id)
        sequences_mapping[peptide_precursor.sequence] = tmp

    for transition in targexp.getTransitions():

        # Fill the precursor mapping
        tmp = precursors_mapping.get(transition.getPeptideRef(), [])
        if transition.getPeptideRef() not in tmp:
            tmp.append(transition.getNativeID())
        precursors_mapping[transition.getPeptideRef()] = tmp

def inferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping,
                 sequences_mapping, protein_mapping, verbose=False, throwOnMismatch=False, fileType=None):
        
    """ Infers a mapping between raw chromatogram files (mzML) and processed feature TSV files

    Usually on feature file can contain multiple aligned runs and maps to
    multiple chromatogram files (mzML). This function will try to guess the
    original name of the mzML based on the align_origfilename column in the
    TSV. Note that both files have some typical endings that are _not_ shared,
    these are generally removed before comparison.

    Only an excact match is allowed.
    """
    import csv, os

    if fileType == "simple":
        return simpleInferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping, sequences_mapping, protein_mapping)
    elif fileType == "traml":
        return tramlInferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping, sequences_mapping, protein_mapping)

    for file_nr, f in enumerate(aligned_pg_files):
        header_dict = {}
        if f.endswith('.gz'):
            import gzip 
            filehandler = gzip.open(f,'rb')
        else:
            filehandler = open(f)
        reader = csv.reader(filehandler, delimiter="\t")
        header = reader.next()
        for i,n in enumerate(header):
            header_dict[n] = i
        if not header_dict.has_key("align_origfilename") or not header_dict.has_key("align_runid"):
            print header_dict
            raise Exception("need column header align_origfilename and align_runid")

        for this_row in reader:

            # Get the mapping ... 
            if header_dict.has_key("FullPeptideName") and \
              header_dict.has_key("Charge") and \
              header_dict.has_key("aggr_Fragment_Annotation"):
                transitions = this_row[ header_dict["aggr_Fragment_Annotation"] ].split(";")
                if len(transitions[-1]) == 0:
                    transitions = transitions[:-1]
                peptide_name = this_row[header_dict["FullPeptideName"]]
                charge_state = this_row[header_dict["Charge"]]
                key = peptide_name + "/" + charge_state
                precursors_mapping [ key ] = transitions
                mapped_precursors = sequences_mapping.get( peptide_name, [] )
                mapped_precursors.append(key)
                sequences_mapping[peptide_name] = [ key ]

                if "ProteinName" in header_dict:
                    protein_name = this_row[header_dict["ProteinName"]]

                    tmp = protein_mapping.get(protein_name, [])
                    if peptide_name not in tmp:
                        tmp.append(peptide_name)
                    protein_mapping[protein_name] = tmp



            # 1. Get the original filename (find a non-NA entry) and the corresponding run id
            if len(this_row) == 0: 
                continue
            if this_row[ header_dict["align_origfilename"] ] == "NA":
                continue
            aligned_id = os.path.basename(this_row[ header_dict["align_runid"] ])
            if aligned_id in mapping:
                continue 


            aligned_fname = os.path.basename(this_row[ header_dict["align_origfilename"] ])

            # 2. Go through all chromatogram input files and try to find
            # one that matches the one from align_origfilename
            for rfile in rawdata_files:

                # 2.1 remove common file endings from the raw data
                rfile_base = os.path.basename(rfile)
                for ending in [".mzML", ".chrom"]:
                    rfile_base = rfile_base.split(ending)[0]

                # 2.2 remove common file endings from the tsv data
                for ending in [".tsv", ".csv", ".xls", "_with_dscore", "_all_peakgroups"]:
                    aligned_fname = aligned_fname.split(ending)[0]

                # 2.3 Check if we have a match
                if aligned_fname == rfile_base:
                    if verbose: 
                        print "- Found match:", os.path.basename(rfile), "->", os.path.basename(this_row[ header_dict["align_origfilename"] ])
                    mapping[aligned_id] = [rfile]

            if not aligned_id in mapping:
                if verbose or throwOnMismatch:
                    print "- No match found for :", aligned_fname, "in any of", \
                            [os.path.basename(rfile) for rfile in rawdata_files]
                    print "- This is generally a very bad sign and you might have " +\
                            "to either rename your files to have matching filenames " +\
                            "or provide an input yaml file describing the matching in detail"
                if throwOnMismatch:
                    raise Exception("Mismatch, alignemnt filename could not be matched to input chromatogram")


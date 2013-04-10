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

#
# The data structures of the reader
#

class MinimalPeakGroup():
    """
    A single peakgroup that is defined by a retention time in a chromatogram
    of multiple transitions. Additionally it has an fdr_score and it has an
    aligned RT (e.g. retention time in normalized space).
    A peakgroup can be selected for quantification or not.
    
    Note that for performance reasons, the peakgroups are created on-the-fly
    and not stored as objects but rather as tuples in "Peptide".

    Each peak group has a unique id, a score (fdr score usually), a retention
    time as well as a back-reference to the precursor that generated the
    peakgroup.
    In this case, the peak group can also be selected or unselected.
    
    """
    def __init__(self, unique_id, fdr_score, assay_rt, selected, peptide, intensity=None):
      self.id = unique_id
      self.fdr_score = fdr_score
      self.diff_from_assay_seconds = assay_rt 
      self.selected_ = selected
      self.peptide = peptide
      self.intensity_ = intensity
  
    ## Print
    def print_out(self):
        # return self.run.get_id() + "/" + self.get_id() + " " + str(self.get_fdr_score()) + " " + str(self.get_normalized_retentiontime()) + " " + str(self.get_value("RT")) + " " + str(self.get_value("rt_score")) # rt_score = delta iRT
        return self.peptide.run.get_id() + "/" + self.get_feature_id() + " " + str(self.get_fdr_score()) + " " + str(self.get_normalized_retentiontime()) # + " " + str(self.get_value("RT")) + " " + str(self.get_value("rt_score")) # rt_score = delta iRT

    ## Getters 
    def get_fdr_score(self):
      return self.fdr_score

    def get_feature_id(self):
      return self.id

    def get_normalized_retentiontime(self):
      return self.diff_from_assay_seconds

    def get_intensity(self):
      return self.intensity_

    ## Select / De-select peakgroup
    def is_selected(self):
        return self.selected_

    ## Setters
    def select_this_peakgroup(self):
        self.selected_ = True
        self.peptide.select_pg(self.id)

    def unselect_this_peakgroup(self):
        self.selected_ = False
        self.peptide.unselect_pg(self.id)

class Precursor():
    """
    A collection of peakgroups that belong to the same precursor and belong
    to one run.
    A peptide can return its best transition group, the selected peakgroup,
    or can return the transition group that is closest to a given iRT time.
    Its id is the transition_group_id (e.g. the id of the chromatogram)
    
    For memory reasons, we store all information about the peakgroup in a
    tuple (invariable). This tuple contains a unique feature id, a score and
    a retention time. Additionally, we also store, whether the feature was
    selected or not.

    A peakgroup has the following attributes: 
        - an identifier that is unique among all other precursors 
        - a set of peakgroups 
        - a backreference to the run it belongs to
    """
    def __init__(self, this_id, run):
        self.id = this_id  
        self.peakgroups = []
        self.run = run
        self.peakgroups_ = []
        self.selected_ = []
        self._decoy = False
  
    def __str__(self):
        return "%s (run %s)" % (self.id, self.run)

    def add_peakgroup_tpl(self, pg_tuple, tpl_id):
        # peakgroup_tuple = (thisid, fdr_score, diff_from_assay_seconds, intensity)
        assert self.id == tpl_id # Check that the peak group is added to the correct precursor
        assert len(pg_tuple) == 4
        self.peakgroups_.append(pg_tuple)
        self.selected_.append(False)

    def get_id(self):
        return self.id 
  
    def get_run_id(self):
      return self.run.get_id()
    
    def get_decoy(self):
        return self._decoy

    def set_decoy(self, decoy):
        if decoy == "FALSE":
            self._decoy = False
        elif decoy == "TRUE":
            self._decoy = True
        else:
            raise Exception("Unknown decoy classifier '%s', please check your input data!" % decoy)
  
    # store information about the peakgroup - tuples (e.g. whether they are selected)
    def select_pg(self, this_id):
        pg_id = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == this_id]
        assert len(pg_id) == 1
        self.selected_[pg_id[0]] = True

    def unselect_pg(self, id):
        pg_id = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == this_id]
        assert len(pg_id) == 1
        self.selected_[pg_id[0]] = False

    def get_best_peakgroup(self):
        if len(self.peakgroups_) == 0: return None
        best_score = self.peakgroups_[0][1]
        result = self.peakgroups_[0]
        for peakgroup in self.peakgroups_:
            if peakgroup[1] <= best_score:
                best_score = peakgroup[1]
                result = peakgroup
        index = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == result[0]][0]
        return MinimalPeakGroup(result[0], result[1], result[2], self.selected_[index], self, result[3])

    def get_selected_peakgroup(self):
      # return the selected peakgroup of this peptide, we can only select 1 or
      # zero groups per chromatogram!
      selected = [i for i,pg in enumerate(self.selected_) if pg]
      assert len(selected) < 2
      if len(selected) == 1:
        index = selected[0]
        result = self.peakgroups_[index]
        return MinimalPeakGroup(result[0], result[1], result[2], self.selected_[index], self, result[3])
      else: 
          return None

    def get_all_peakgroups(self):
        for index, result in enumerate(self.peakgroups_):
            yield MinimalPeakGroup(result[0], result[1], result[2], self.selected_[index], self, result[3])
  
    def find_closest_in_iRT(self, delta_assay_rt):
      result = min(self.peakgroups_, key=lambda x: abs(float(x[2]) - float(delta_assay_rt)))
      index = [i for i,pg in enumerate(self.peakgroups_) if pg[0] == result[0]][0]
      return MinimalPeakGroup(result[0], result[1], result[2], self.selected_[index], self, result[3])

class Run():
    """
    One single SWATH run that contains peptides (chromatograms) 
    It has a unique id and stores the headers from the csv

    A run has the following attributes: 
        - an identifier that is unique to this run
        - a filename where it originally came from
        - a dictionary of precursors, accessible through a dictionary
    """

    def __init__(self, header, header_dict, runid, orig_filename):
        self.header = header
        self.header_dict = header_dict
        self.runid = runid
        self.orig_filename = orig_filename
        self.all_peptides = {}
  
    def get_id(self):
        return self.runid
  
    def get_best_peaks(self):
        result = []
        for k, peptide in self.all_peptides.iteritems():
          result.append(peptide.get_best_peakgroup())
        return result
  
    def get_best_peaks_with_cutoff(self, cutoff):
        return [p for p in self.get_best_peaks() if p.get_fdr_score() < cutoff]
  
    def get_all_trgroups(self, cutoff):
      above_cutoff = []
      for k,peak in self.all_peptides.iteritems():
        if peak.get_fdr_score() < cutoff:
            above_cutoff.append(peak)
      return above_cutoff
  
    def get_peptide(self, id):
        try:
          return self.all_peptides[id]
        except KeyError:
          # this run has no peakgroup for that peptide
          return None

    def __iter__(self):
        for peptide in self.all_peptides.values():
            yield peptide


#
# The Readers of the Scoring files
#

class SWATHScoringReader:

    def __init__(self):
        raise Exception("Abstract class")

    def parse_row(self, run, this_row, do_realignment):
        raise Exception("Abstract method")

    @staticmethod
    def newReader(infiles, filetype):
        """Factory to create a new reader"""
        if filetype  == "openswath": 
            return OpenSWATH_SWATHScoringReader(infiles)
        elif filetype  == "peakview": 
            return Peakview_SWATHScoringReader(infiles)
        else:
            raise Exception("Unknown filetype '%s', allowed types are %s" % (decoy, str(filetypes) ) )

    def parse_files(self, do_realignment):
      print "Parsing input files"

      from sys import stdout
      import csv
      runs = []
      for file_nr, f in enumerate(self.infiles):
        stdout.write("\rReading %s" % str(f))
        stdout.flush()
        header_dict = {}
        filehandler = open(f)
        reader = csv.reader(filehandler, delimiter="\t")
        header = reader.next()
        for i,n in enumerate(header):
          header_dict[n] = i
        stdout.write("\rReading file %s" % (str(f)) )
        stdout.flush()
        # There may be multiple runs in one csv file, we use the run number as
        # well as the file number to find unique runs
        for this_row in reader:
            # check if we have a new run
            runnr = this_row[header_dict[self.run_id_name]]
            runid = runnr + "_" + str(file_nr)
            current_run = [r for r in runs if r.get_id() == runid]
            if len(current_run) == 0:
                current_run = Run(header, header_dict, runid, f)
                runs.append(current_run)
            else: 
                assert len(current_run) == 1
                current_run = current_run[0]
            # Unfortunately, since we are using csv tell() will not work...
            # print "parse row at", filehandler.tell()
            self.parse_row(current_run, this_row, do_realignment)
      # Here we check that each run indeed has a unique id
      assert len(set([r.get_id() for r in runs])) == len(runs) # each run has a unique id
      stdout.write("\r\r\n") # clean up
      print "Found %s runs" % len(runs)
      return runs

class OpenSWATH_SWATHScoringReader(SWATHScoringReader):

    def __init__(self, infiles):
        self.infiles = infiles
        self.run_id_name = "run_id"

    def parse_row(self, run, this_row, do_realignment):
        decoy_name = "decoy"
        fdr_score_name = "m_score"
        unique_peakgroup_id_name = "transition_group_id"
        diff_from_assay_in_sec_name = "delta_rt"
        run_id_name = "run_id"
        protein_id_col = "ProteinName"
        sequence_col = "Sequence"
        unique_feature_id_name = "id"
        intensity_name = "Intensity"
        decoy = "FALSE"

        # use the aligned retention time if it is available!
        if "aligned_rt" in run.header_dict: 
            diff_from_assay_in_sec_name = "aligned_rt" ## use this if it is present
        # if we want to re-do the re-alignment, we just use the "regular" retention time
        if do_realignment: 
            diff_from_assay_in_sec_name = "RT"

        trgr_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        protein_name = this_row[run.header_dict[protein_id_col]]
        sequence = this_row[run.header_dict[sequence_col]]
        thisid = this_row[run.header_dict[unique_feature_id_name]]
        fdr_score = float(this_row[run.header_dict[fdr_score_name]])
        diff_from_assay_seconds = float(this_row[run.header_dict[diff_from_assay_in_sec_name]])
        unique_peakgroup_id = this_row[run.header_dict[unique_peakgroup_id_name]]
        intensity = float(this_row[run.header_dict[intensity_name]])
        if "decoy" in run.header_dict:
            decoy = this_row[run.header_dict[decoy_name]]
        run_id = int(this_row[run.header_dict[run_id_name]])

        if not run.all_peptides.has_key(trgr_id):
          p = Precursor(trgr_id, run)
          p.protein_name = protein_name
          p.sequence = sequence
          p.run_id = run_id
          p.set_decoy(decoy)
          run.all_peptides[trgr_id] = p
        peakgroup_tuple = (thisid, fdr_score, diff_from_assay_seconds, intensity)
        run.all_peptides[trgr_id].add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id)

class Peakview_SWATHScoringReader(SWATHScoringReader):

    def __init__(self, infiles):
        self.infiles = infiles
        self.run_id_name = "Sample"

    def parse_row(self, run, this_row, do_realignment):
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
        if do_realignment: 
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

        if not run.all_peptides.has_key(trgr_id):
          p = Precursor(trgr_id, run)
          p.protein_name = protein_name
          p.sequence = sequence
          p.run_id = run_id
          p.set_decoy(decoy)
          run.all_peptides[trgr_id] = p
          if verb: print "add peptide", trgr_id
        peakgroup_tuple = (thisid, fdr_score, diff_from_assay_seconds,intensity)
        if verb: print "append tuple", peakgroup_tuple
        run.all_peptides[trgr_id].add_peakgroup_tpl(peakgroup_tuple, unique_peakgroup_id)



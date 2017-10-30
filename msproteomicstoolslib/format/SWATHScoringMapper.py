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
import os

maxInt = maxsize
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

"""
Doc :
    A class to read and map the SWATH from Peakview and OpenSWATH
"""


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
        print("\nError!")
        print("Could not import pyOpenMS while trying to load a TraML file - please make sure pyOpenMS is installed.")
        print("pyOpenMS is available from https://pypi.python.org/pypi/pyopenms")
        print()
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

def mapRow(this_row, header_dict, precursors_mapping, sequences_mapping, protein_mapping):
    # Get the mapping ... 
    if "FullPeptideName" in header_dict and \
      "Charge" in header_dict and \
      "aggr_Fragment_Annotation" in header_dict:
        transitions = this_row[ header_dict["aggr_Fragment_Annotation"] ].split(";")
        pr_transitions = []
        if "aggr_prec_Fragment_Annotation" in header_dict:
            pr_transitions = this_row[ header_dict["aggr_prec_Fragment_Annotation"] ].split(";")
        if len(transitions[-1]) == 0:
            transitions = transitions[:-1]
        if len(pr_transitions) > 0 and len(pr_transitions[-1]) == 0:
            pr_transitions = pr_transitions[:-1]
        peptide_name = this_row[header_dict["FullPeptideName"]]
        charge_state = this_row[header_dict["Charge"]]
        key = peptide_name + "/" + charge_state
        prkey = peptide_name + "/" + charge_state + "_pr"
        precursors_mapping [ key ] = transitions
        precursors_mapping [ prkey ] = pr_transitions
        mapped_precursors = sequences_mapping.get( peptide_name, [] )
        mapped_precursors.extend([key, prkey])
        sequences_mapping[peptide_name] = mapped_precursors #  = [ key, prkey ]

        if "ProteinName" in header_dict:
            protein_name = this_row[header_dict["ProteinName"]]

            tmp = protein_mapping.get(protein_name, [])
            if peptide_name not in tmp:
                tmp.append(peptide_name)
            protein_mapping[protein_name] = tmp

def getAlignedFilename(this_row, header_dict):

    if this_row[ header_dict["align_origfilename"] ] == "NA":
        return None, None
    aligned_id = os.path.basename(this_row[ header_dict["align_runid"] ])

    aligned_fname = os.path.basename(this_row[ header_dict["align_origfilename"] ])

    # allow cross-OS transfer of files
    tmp = aligned_fname.split("\\")
    aligned_fname = tmp[-1]
    tmp = aligned_fname.split("/")
    aligned_fname = tmp[-1]
    return aligned_fname, aligned_id

def sqlInferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping,
                 sequences_mapping, protein_mapping, verbose=False, throwOnMismatch=False, fileType=None):

    import sqlite3, csv, os

    sqlfile_map = []
    for k, filename in enumerate(rawdata_files):
        conn = sqlite3.connect(filename)
        c = conn.cursor()
        d = list(c.execute("SELECT ID, FILENAME, NATIVE_ID FROM RUN"))
        assert len(d) == 1
        filename = d[0][1]
        sqlfile_map.append([k, 0, os.path.basename(filename)])
        mapping[k] = [None]

    nomatch_found = set([])
    for file_nr, f in enumerate(aligned_pg_files):
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
        if not "align_origfilename" in header_dict or not "align_runid" in header_dict:
            print (header_dict)
            raise Exception("need column header align_origfilename and align_runid")

        for this_row in reader:

            # Get the transition mapping ... 
            mapRow(this_row, header_dict, precursors_mapping, sequences_mapping, protein_mapping)

            # 1. Get the original filename (find a non-NA entry) and the corresponding run id
            aligned_fname, aligned_id = getAlignedFilename(this_row, header_dict)

            # 2. Go through all sql input files and try to find
            # one that matches the one from align_origfilename
            for sqlfile, runid, rfile in sqlfile_map:

                # 2.1 remove common file endings from the raw data
                rfile_base = rfile
                for ending in [".gz", ".mzML", ".chrom", ".sqMass"]:
                    rfile_base = rfile_base.split(ending)[0]

                # 2.2 remove common file endings from the tsv data
                for ending in [".tsv", ".csv", ".xls", "_with_dscore", "_all_peakgroups"]:
                    aligned_fname = aligned_fname.split(ending)[0]

                # 2.3 Check if we have a match
                # print aligned_fname, rfile_base
                if aligned_fname == rfile_base:
                    if False:
                        print("- Found match:", os.path.basename(rfile), "->", os.path.basename(this_row[ header_dict["align_origfilename"] ]))
                    mapping[sqlfile][runid] = aligned_id

def inferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping,
                 sequences_mapping, protein_mapping, verbose=False, throwOnMismatch=False, fileType=None):
        
    """ Infers a mapping between raw chromatogram files (mzML) and processed feature TSV files

    Usually one feature file can contain multiple aligned runs and maps to
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
    elif fileType == "sqmass":
        return sqlInferMapping(rawdata_files, aligned_pg_files, mapping, precursors_mapping, sequences_mapping, protein_mapping)

    nomatch_found = set([])
    for file_nr, f in enumerate(aligned_pg_files):
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
        if not "align_origfilename" in header_dict or not "align_runid" in header_dict:
            print (header_dict)
            raise Exception("need column header align_origfilename and align_runid")

        for this_row in reader:

            if len(this_row) == 0: 
                continue

            # Get the transition mapping ... 
            mapRow(this_row, header_dict, precursors_mapping, sequences_mapping, protein_mapping)

            # 1. Get the original filename (find a non-NA entry) and the corresponding run id
            aligned_fname, aligned_id = getAlignedFilename(this_row, header_dict)
            if aligned_id is None or aligned_id in mapping:
                continue 

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
                        print("- Found match:", os.path.basename(rfile), "->", os.path.basename(this_row[ header_dict["align_origfilename"] ]))
                    mapping[aligned_id] = [rfile]

            if not aligned_id in mapping:
                if verbose or throwOnMismatch:
                    nomatch_found.update( [aligned_fname] )
                if throwOnMismatch:
                    raise Exception("Mismatch, alignment filename could not be matched to input chromatogram")

        if verbose:
            print("- No match found for :", list(nomatch_found), "in any of", \
              [os.path.basename(rfile) for rfile in rawdata_files])
            print("- This may be a bad sign if you expected a match here. You might have " +\
                    "to either rename your files to have matching filenames " +\
                    "or provide an input yaml file describing the matching in detail.")


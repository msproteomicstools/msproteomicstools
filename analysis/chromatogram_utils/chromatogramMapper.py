#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
=========================================================================
        DIAlignPy -- Alignment of Targeted Mass Spectrometry Runs
=========================================================================

<Shubham Gupta chromatogramMapper.py>
Copyright (C) 2020 Shubham Gupta

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

--------------------------------------------------------------------------
$Maintainer: Shubham Gupta$
$Authors: Shubham Gupta$
--------------------------------------------------------------------------
"""

from pyopenms import OnDiscMSExperiment

class mz_pointer():
    def __init__(self, filename):
        # Establish connection to mzML files
        mz = OnDiscMSExperiment()
        mz.openFile(filename)
        self.filename = filename
        self.pointer = mz
        # Get maximum number of chromatogram
        self.meta_data = mz.getMetaData()
        self.maxIndex = len(self.meta_data.getChromatograms()) - 1 

    def extractXIC_group_(self, chromIndices):
        """ Extract chromatograms for chromatrogram indices """
        XIC_group = [None]*len(chromIndices)
        for i in range(len(chromIndices)):
            # Chromatogram is a tuple (time_array, intensity_array) of numpy array
            # mz.getChromatogramById(chromIndices[i]).getIntensityArray()
            # mz.getChromatogramById(chromIndices[i]).getTimeArray()
            if chromIndices[i] > self.maxIndex:
                print("Warning : Invalid index {} in run {}. Skipping!".format(chromIndices[i], self.filename))
                continue
            XIC_group[i] = self.pointer.getChromatogram(chromIndices[i]).get_peaks()
        return XIC_group

class mzml_accessors():
    """ Establish connection to mzML files """
    def __init__(self, runs, MStoFeature):
        # TODO Should I add the mz accessor to the Run object ?
        chrom_file_accessor = {}
        for run in runs:
            mzml_file = MStoFeature[run.get_openswath_filename()][0]
            chrom_file_accessor[run.get_id()] = mz_pointer(mzml_file)
        self.pointers = chrom_file_accessor
        self.runs = runs
        self.invalid_chromIndex = -1
        self.run_chromIndex_map = {}
    
    def extractXIC_group(self, run, prec_id):
        chrom_ids = self.run_chromIndex_map[run.get_id()][prec_id]
        if self.invalid_chromIndex in chrom_ids:
            print("Can't extract XICs for precursor {} from run {}. Skipping!".format(prec_id, run.get_id()))
            return None
        else:
            return self.pointers[run.get_id()].extractXIC_group_(chrom_ids)

    def get_precursor_to_chromID(self, precursor_to_transitionID, invalid_chromIndex = -1):
        self.invalid_chromIndex = invalid_chromIndex
        # Get precursor to chromatogram indices
        for run in self.runs:
            meta_data = self.pointers[run.get_id()].meta_data
            precursor_to_chromIndex = {}
            # Initialize chromatogram indices with -1
            for prec in precursor_to_transitionID:
                precursor_to_chromIndex[prec] = [invalid_chromIndex]*len(precursor_to_transitionID[prec][1])
            # Iterate through Native IDs in mzML file
            for i in range(meta_data.getNrChromatograms()):
                nativeID = int(meta_data.getChromatogram(i).getNativeID())
                # Check if Native ID matches to any of our transition IDs.
                for prec in precursor_to_transitionID:
                    if nativeID in precursor_to_transitionID[prec][1]:
                        # Find the index of transition ID and insert chromatogram index
                        index = precursor_to_transitionID[prec][1].index(nativeID)
                        precursor_to_chromIndex[prec][index] = i
            self.run_chromIndex_map[run.get_id()] = precursor_to_chromIndex

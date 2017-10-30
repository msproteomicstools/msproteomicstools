#!/usr/bin/python
# -*- coding: utf-8 -*-
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

CHROMTYPES = {
    0 : "Protein", 
    1 : "Peptide", 
    2 : "Precursor", 
    3 : "Transition"
} 

CHROMTYPES_r = dict([ (v,k) for k,v in CHROMTYPES.iteritems()])

class ChromatogramTransition(object):
    """
    Internal tree structure object representing one row in the in the left side tree.

    This is the bridge between the view and the data model

    Pointers to objects of :class:`.ChromatogramTransition` are passed to
    callback functions when the selection of the left side tree changes. The
    object needs to have store information about all the column present in the
    rows (PeptideSequence, Charge, Name) which are requested by the
    :class:`.PeptideTree` model.

    Also it needs to know how to access the raw data as well as meta-data for a
    certain transition.  This is done through getData, getLabel etc.
    """

    def __init__(self, name, charge, subelements, peptideSequence=None, fullName=None, datatype="Precursor"):
        self._name = name
        self._charge = charge
        self._fullName = fullName
        self._peptideSequence = peptideSequence
        self._subelements = subelements
        self.mytype = CHROMTYPES_r[datatype]

    def getSubelements(self):
        return self._subelements

    def getPeptideSequence(self):
        if self._peptideSequence is None:
            return self.getName()
        return self._peptideSequence

    def getName(self):
        """
        Get name of precursor

        Returns
        -------
        str:
            Name of precursor
        """
        return self._name

    def getCharge(self):
        """
        Get charge of precursor

        Returns
        -------
        int:
            Charge
        """
        return self._charge

    def getType(self):
        return CHROMTYPES[self.mytype]

    def getData(self, run):
        """
        Get raw data for a certain object

        If we have a single precursors or a peptide with only one precursor, we
        show the same data as for the precursor itself. For a peptide with
        multiple precursors, we show all precursors as individual curves. For a
        single transition, we simply plot that transition.

        Parameters
        ----------
        run : :class:`.SwathRun` or :class:`.SqlSwathRun`
            SwathRun object which will be used to retrieve data

        Returns
        -------
        list of pairs (timearray, intensityarray): 
            Returns the raw data of the chromatograms for a given run. The
            dataformat is a list of transitions and each transition is a pair
            of (timearray,intensityarray)
        """

        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_data_for_precursor(self.getName()) 

        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_data_for_precursor(prec[0]) 
            else:

                # Peptide view with multiple precursors
                # -> Sum up the data for all individual precursors
                final_data = []
                for p in prec:
                    timedata = None
                    intdata = None
                    import numpy
                    for data in run.get_data_for_precursor(p):
                        if timedata is None:
                            timedata = numpy.array(data[0])
                            intdata = numpy.array(data[1])
                        else:
                            intdata = intdata + numpy.array(data[1])
                    final_data.append( [timedata, intdata] )

                return final_data

        elif CHROMTYPES[self.mytype] == "Transition" :
            return run.get_data_for_transition(self.getName()) 

        return [ [ [0], [0] ] ]

    def getRange(self, run):
        """
        Get the data range (leftWidth/rightWidh) for a specific run

        Parameters
        ----------
        run : :class:`.SwathRun`
            SwathRun object which will be used to retrieve data

        Returns
        -------
        list of float:
            A pair of floats representing the data range (leftWidth/rightWidh) for a specific run
        """

        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_range_data(self.getName()) 

        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_range_data(prec[0]) 

        elif CHROMTYPES[self.mytype] == "Transition" :
            # TODO
            return [ [0,0] ]

        return [ [0,0] ]

    def getProbScore(self, run):
        """
        Get the probabilistic score for a specific run and current precursor

        Parameters
        ----------
        run : :class:`.SwathRun`
            SwathRun object which will be used to retrieve data

        Returns
        -------
        float:
            The probabilistic score for a specific run and current precursor
        """

        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_score_data(self.getName()) 

        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_score_data(prec[0]) 
            else: 
                # For multiple precursors, the probability score is not defined 
                return None

        elif CHROMTYPES[self.mytype] == "Transition" :
            return None

        return None

    def getIntensity(self, run):
        """
        Get the intensity for a specific run and current precursor

        Parameters
        ----------
        run : :class:`.SwathRun`
            SwathRun object which will be used to retrieve data

        Returns
        -------
        float:
            The intensity for a specific run and current precursor
        """

        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_intensity_data(self.getName()) 

        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_intensity_data(prec[0]) 
            else: 
                # For multiple precursors, the intensity is currently not computed 
                return None

        elif CHROMTYPES[self.mytype] == "Transition" :
            return None

        return None

    def getAssayRT(self, run):
        """
        Get the intensity for a specific run and current precursor

        Parameters
        ----------
        run : :class:`.SwathRun`
            SwathRun object which will be used to retrieve data

        Returns
        -------
        float:
            The intensity for a specific run and current precursor
        """

        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_assay_data(self.getName()) 

        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_assay_data(prec[0]) 
            else: 
                # For multiple precursors, the intensity is currently not computed 
                return None

        elif CHROMTYPES[self.mytype] == "Transition" :
            return None

        return None

    def getLabel(self, run):
        """
        Get the labels for a curve (corresponding to the raw data from getData
        call) for a certain object.

        If we have a single precursors or a peptide with only one precursor, we
        show the same data as for the precursor itself. For a peptide with
        multiple precusors, we show all precursors as individual curves. For a
        single transition, we simply plot that transition.

        Parameters
        ----------
        run : :class:`.SwathRun`
            SwathRun object which will be used to retrieve data

        Returns
        -------
        list of str:
            The labels to display for each line in the graph
        """

        if CHROMTYPES[self.mytype] == "Precursor" :
            return run.get_transitions_for_precursor_display(self.getName())

        elif CHROMTYPES[self.mytype] == "Peptide" :
            prec = run.get_precursors_for_sequence(self.getName())
            if len(prec) == 1:
                return run.get_transitions_for_precursor_display(prec[0])
            else:
                # Peptide view with multiple precursors
                return prec

        elif CHROMTYPES[self.mytype] == "Transition" :
            return [self.getName()]

        return [ "" ]


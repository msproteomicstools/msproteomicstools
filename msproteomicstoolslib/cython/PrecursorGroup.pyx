# distutils: language = c++
# cython: c_string_encoding=ascii  # for cython>=0.19
# encoding: latin-1
cimport cython
cimport libc.stdlib
cimport numpy as np

cdef class CyPrecursorGroup(object):
    """See :class:`.PrecursorGroup` for a description.

    This implementation is pure Cython.

    Attributes:
        - self.peptide_group_label_: Identifier or precursor group 
        - self.run_: Reference to the :class:`.Run` where this PrecursorGroup is from
        - self.precursors_: List of :class:`.CyPrecursorWrapperOnly`
    """

    __slots__ = ["peptide_group_label_", "run_", "precursors_"]

    cdef libcpp_string peptide_group_label_ 
    cdef object run_
    cdef list precursors_ # list of peptide precursors (CyPrecursor)

    def __init__(self, peptide_group_label, run):
        self.peptide_group_label_ = peptide_group_label  
        self.run_ = run
        self.precursors_ = []

    def __str__(self):
        return "PrecursorGroup %s" % (self.getPeptideGroupLabel())

    ### def __lt__(self, other):

    ###     if self.run_.get_id() == other.run_.get_id():
    ###         return self.getPeptideGroupLabel() > other.getPeptideGroupLabel()
    ###     else:
    ###         return self.run_.get_id() > other.run_.get_id()

    def __iter__(self):
        ## # print ("Cy - iter through pre")
        ## for precursor in self.precursors_:
        ##     print ("Will iterate and provide ", type(precursor))

        # print ("Cy - iter through pre")
        for precursor in self.precursors_:
            ### print ("Cy - test1 ")
            ### print (type(precursor))
            yield precursor
            ## print ("Cy - test2 ")
        # print ("Cy -- done")

    def __classInvariant__(self):
        if len(self.precursors_) > 0:
            # All precursor sequences should all be equal to the first sequence
            assert(all( [precursor.getSequence() == self.precursors_[0].getSequence() for precursor in self.precursors_] )) 
        return True

    # @class_invariant(__classInvariant__)
    def getPeptideGroupLabel(self):
        """
        getPeptideGroupLabel(self)
        Get peptide group label
        """
        return self.peptide_group_label_
  
    # @class_invariant(__classInvariant__)
    def addPrecursor(self, precursor):
        """
        addPrecursor(self, precursor)
        Add precursor to peptide group
        """
        precursor.set_precursor_group( self )
        self.precursors_.append(precursor)

    # @class_invariant(__classInvariant__)
    def getPrecursor(self, curr_id):
        """
        getPrecursor(self, curr_id)
        Get the precursor for the given transition group id
        """
        for precursor in self:
            if precursor.get_id() == curr_id:
                return precursor
        return None

    # @class_invariant(__classInvariant__)
    def getAllPrecursors(self):
        """
        getAllPrecursors(self)
        Return a list of all precursors in this precursor group
        """
        return list(self)

    # @class_invariant(__classInvariant__)
    def getAllPeakgroups(self):
        """
        getAllPeakgroups(self)
        Generator of all peakgroups attached to the precursors in this group
        """
        for pr in self.precursors_:
            for pg in pr.get_all_peakgroups():
                yield pg

    # @class_invariant(__classInvariant__)
    def getOverallBestPeakgroup(self):
        """
        getOverallBestPeakgroup(self)
        Get the best peakgroup (by fdr score) of all precursors contained in this precursor group
        """
        allpg = list(self.getAllPeakgroups())
        if len(allpg) == 0:
            return None

        minscore = min([pg.get_fdr_score() for pg in allpg])
        return [pg for pg in allpg if pg.get_fdr_score() <= minscore][0]

    def get_decoy(self):
        """
        Whether the current peptide is a decoy or not

        Returns:
            decoy(bool): Whether the peptide is decoy or not
        """
        if len(self.precursors_) == 0:
            return False

        return self.precursors_[0].get_decoy()

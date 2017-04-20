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

import MSData

from PyQt4 import QtCore 
from PyQt4.QtCore import Qt, QModelIndex
from TreeModels import TreeNode
from TreeModels import TreeModel

# A TreeNode element
class PeptideTreeNode(TreeNode):
    """
    Implementation of a node in the right-hand peptide tree in the GUI
    """

    def __init__(self, ref, parent, row):
        self.ref = ref
        TreeNode.__init__(self, parent, row)

    def _getChildren(self):
        return [PeptideTreeNode(elem, self, index)
            for index, elem in enumerate(self.ref.getSubelements())]

class PeptideTree(TreeModel):
    """
    Implementation of a tree in the right-hand peptide tree in the GUI
    """

    def __init__(self, rootElements, firstColumnName="Peptide Sequence"):

        self.rootElements = rootElements
        self.first_column_name_ = firstColumnName
        TreeModel.__init__(self)

    def _getRootNodes(self):
        """
        Return root nodes (top-level nodes)
        """
        return [PeptideTreeNode(elem, None, index)
            for index, elem in enumerate(self.rootElements)]

    def columnCount(self, parent):
        """
        Returns how many columns we have
        """
        return 3

    def data(self, index, role):
        """
        Get data for a specific index (and role)

        Currently supported role is only Qt.DisplayRole (for displaying the
        tree). The three columns are:

            - Compound name (generally peptide sequence or compound sum formula)
            - Charge
            - Name

        Parameters
        ----------
        index : QModelIndex 
            Index of the element to be accessed
        role : Qt::ItemDataRole
            Item role to be used (only Qt.DisplayRole supported)
        """
        if not index.isValid():
            return None

        node = index.internalPointer()
        if role == Qt.DisplayRole and index.column() == 0:
            # return QtCore.QVariant( node.ref.getPeptideSequence() )
            return node.ref.getPeptideSequence()
        if role == Qt.DisplayRole and index.column() == 1:
            # return QtCore.QVariant( node.ref.getCharge() )
            return node.ref.getCharge()
        if role == Qt.DisplayRole and index.column() == 2:
            # return QtCore.QVariant( node.ref.getName() )
            return node.ref.getName()

        return None

    def headerData(self, section, orientation, role):
        """
        Get header data (column header) for a specific index (and role)

        The three columns are:
            - Peptide Sequence
            - Charge
            - Name

        Note that the user can set the name of the first column name manually
        in order to accomodate for other data (e.g. metabolomics) where
        "Peptide Sequence" would not make sense.
        """
        if section == 0:
            return self.first_column_name_
        if section == 1:
            return 'Charge'
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 2:
            return 'Full Name'
        return None

    def set_precursor_tree_structure(self, data, sortData = True):
        """
        Initialize tree structure with data from :meth:`.MSData.get_precursor_tree`

        The tree is initialized by giving it a pointer to the root element(s)

        Parameters
        ----------
        data : list of :class:`.ChromatogramTransition`:
            Root element(s) for the peptide tree
        sortData : bool
            Whether to sort data
        """

        if sortData:
            data.sort(lambda x,y: cmp(x.getName(), y.getName() ))

        # first delete all rows
        parent = QModelIndex()
        self.beginRemoveRows(parent, 0, len(self.rootElements) )
        self.rootElements = []
        self.endRemoveRows()

        # now add the new rows
        parent = QModelIndex()
        self.beginInsertRows(parent, 0, len(data) )
        self.rootElements = data
        self.reset()
        self.endInsertRows()

        # initialize super method again
        self.initialize()


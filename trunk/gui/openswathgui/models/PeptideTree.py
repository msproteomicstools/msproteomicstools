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
    def __init__(self, ref, parent, row):
        self.ref = ref
        TreeNode.__init__(self, parent, row)

    def _getChildren(self):
        return [PeptideTreeNode(elem, self, index)
            for index, elem in enumerate(self.ref.getSubelements())]

class PeptideTree(TreeModel):
    def __init__(self, rootElements):
        self.rootElements = rootElements
        TreeModel.__init__(self)

    def _getRootNodes(self):
        return [PeptideTreeNode(elem, None, index)
            for index, elem in enumerate(self.rootElements)]

    def columnCount(self, parent):
        return 3

    def data(self, index, role):
        if not index.isValid():
            return None
        node = index.internalPointer()
        if role == Qt.DisplayRole and index.column() == 0:
            return QtCore.QVariant( node.ref.getPeptideSequence() )
        if role == Qt.DisplayRole and index.column() == 1:
            return QtCore.QVariant( node.ref.getCharge() )
        if role == Qt.DisplayRole and index.column() == 2:
            return QtCore.QVariant( node.ref.getName() )
        return None

    def headerData(self, section, orientation, role):
        if section == 0:
            return 'Peptide Sequence'
        if section == 1:
            return 'Charge'
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 2:
            return 'Full Name'
        return None

    def set_precursor_tree_structure(self, data):

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


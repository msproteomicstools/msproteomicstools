#!/usr/bin/python
# -*- coding: utf-8 -*-

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
            for index, elem in enumerate(self.ref.subelements)]

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
            return QtCore.QVariant( node.ref.charge )
        if role == Qt.DisplayRole and index.column() == 2:
            return QtCore.QVariant( node.ref.name )
        return None

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 0:
            return 'Peptide Sequence'
        if orientation == Qt.Horizontal and role == Qt.DisplayRole \
            and section == 1:
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
        self.endInsertRows()

        # initialize super method again
        self.initialize()

    def setHorizontalHeaderLabels(self, parent):
        pass

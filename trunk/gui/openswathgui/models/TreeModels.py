#!/usr/bin/python
# -*- coding: utf-8 -*-

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt, QModelIndex

## Generic tree models
# from http://www.hardcoded.net/articles/using_qtreeview_with_qabstractitemmodel.htm
class TreeNode(object):
    """
    Generic model of a tree node

    Adopted from http://www.hardcoded.net/articles/using_qtreeview_with_qabstractitemmodel.htm

    See :class:`.PeptideTreeNode` for implementation.
    """
    def __init__(self, parent, row):
        self.parent = parent
        self.row = row
        self.subnodes = self._getChildren()

    def _getChildren(self):
        """
        Get children of current node
        """
        raise NotImplementedError()


class TreeModel(QtCore.QAbstractItemModel):
    """
    Generic tree model

    Adopted from http://www.hardcoded.net/articles/using_qtreeview_with_qabstractitemmodel.htm

    See parent class http://qt-project.org/doc/qt-5/QAbstractItemModel.html

    See :class:`.PeptideTree` for implementation.
    """
    def __init__(self):
        QtCore.QAbstractItemModel.__init__(self)
        self.rootNodes = self._getRootNodes()

    def initialize(self):
        """
        Initialize tree
        """
        self.rootNodes = self._getRootNodes()

    def _getRootNodes(self):
        raise NotImplementedError()

    def index(self, row, column, parent):
        if not parent.isValid():
            return self.createIndex(row, column, self.rootNodes[row])
        parentNode = parent.internalPointer()
        return self.createIndex(row, column, parentNode.subnodes[row])

    def parent(self, index):
        if not index.isValid():
            return QModelIndex()
        node = index.internalPointer()
        if node.parent is None:
            return QModelIndex()
        else:
            return self.createIndex(node.parent.row, 0, node.parent)

    def reset(self):
        self.rootNodes = self._getRootNodes()
        QtCore.QAbstractItemModel.reset(self)

    def rowCount(self, parent):
        if not parent.isValid():
            return len(self.rootNodes)
        node = parent.internalPointer()
        return len(node.subnodes)


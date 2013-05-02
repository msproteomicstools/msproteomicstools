#!/usr/bin/python
# -*- coding: utf-8 -*-

from PyQt4 import QtGui, QtCore

class PeptidesTreeView( QtGui.QTreeView ):

    def __init__(self):
        super(PeptidesTreeView, self).__init__()

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        # self.customContextMenuRequested.connect(self.openMenu)
        
    def expandMultiElementItems(self):
        """
        Expand all top elements in the tree that have more than one child
        """
        for model_idx in self.iterTopLevelElements(0):
            if len(model_idx.internalPointer().subnodes) > 1:
                self.setExpanded(model_idx, True)

    def iterTopLevelElements(self, column):
        # find (fake) root and model
        root = self.rootIndex()
        m = self.model()
        for i in range(m.rowCount(root)):
            yield m.index(i,column,root)


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


from PyQt4 import QtGui, QtCore

class PeptidesTreeView( QtGui.QTreeView ):
    """
    The Peptide Tree View widget is the view implementation for the left side tree view
    """

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

    def iterAllLevelElements(self, column):
        root = self.rootIndex()
        m = self.model()
        for elem in self._iterAllLevelElements_rec(column, m, root):
            yield elem

    def _iterAllLevelElements_rec(self, column, m, root, level=1):
        for i in range(m.rowCount(root)):
            yield m.index(i,column,root)
            for elem in self._iterAllLevelElements_rec(column, m, m.index(i,column,root), level+1):
                yield elem

    def selectAndScrollTo(self, model_idx):
        if model_idx is None:
            return

        selectionModel = self.selectionModel()
        selectionModel.clearSelection()
        selectionModel.select(model_idx, QtGui.QItemSelectionModel.Select)
        self.setSelectionModel(selectionModel)
        # Now scroll to the item
        self.scrollTo(model_idx, QtGui.QAbstractItemView.PositionAtCenter)


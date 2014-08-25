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

import unittest
import os
from nose.plugins.attrib import attr

import msproteomicstoolslib.algorithms.graphs.graphs as graphs

class TestUnitGraphs(unittest.TestCase):

    def setUp(self):
        # A linear graph
        self.mst1 = [('0_2', '0_3'), ('0_2', '0_0'), ('0_0', '0_1')]

        # A non linear graph
        self.mst3 = [('0_0', '0_1'), ('0_1', '0_2'), ('0_2', '0_3'), ('0_1', '0_4'), ('0_4', '0_5')]
        # an extra node from 0 to 3, this is not an MST
        self.graph3 = [('0_0', '0_1'), ('0_1', '0_2'), ('0_2', '0_3'), ('0_1', '0_4'), ('0_4', '0_5'), ('0_0', '0_3')]

    def test_getAdj_1(self):
        """Test the generation of the adjacency list"""
        adj = graphs.getAdjacencyList(self.mst3)
        expected = {
            '0_4': ['0_1', '0_5'],
            '0_5': ['0_4'], 
            '0_2': ['0_1', '0_3'],
            '0_3': ['0_2'],
            '0_0': ['0_1'], 
            '0_1': ['0_0', '0_2', '0_4']
        }

        for keys in expected:
            self.assertEqual( sorted(expected[keys]), sorted(adj[keys]) )

    def test_getAdj_2(self):
        """Test the generation of the adjacency list"""
        adj = graphs.getAdjacencyList(self.graph3)
        expected = {
            '0_4': ['0_1', '0_5'],
            '0_5': ['0_4'], 
            '0_2': ['0_1', '0_3'],
            '0_3': ['0_2', '0_0'],
            '0_0': ['0_1', '0_3'], 
            '0_1': ['0_0', '0_2', '0_4']
        }

        for keys in expected:
            self.assertEqual( sorted(expected[keys]), sorted(adj[keys]) )

    def test_doBFS_1(self):
        """Test the BFS algorithm"""
        res = list(graphs.doBFS(self.mst3, "0_3"))
        # Note that 0_0 should be before 0_4 in a sorted BFS
        expected = ['0_3', '0_2', '0_1', '0_0', '0_4', '0_5']
        self.assertEqual(expected, res)

    def test_doBFS_2(self):
        """Test the BFS algorithm"""
        res = list(graphs.doBFS(self.graph3, "0_3"))
        # Note that 0_0 should be before 0_2 in a sorted BFS
        expected = ['0_3', '0_0', '0_2','0_1', '0_4', '0_5']
        self.assertEqual(expected, res)

    def test_findShortestMSTPath(self):
        """Test the findShortestMSTPath algorithm"""

        # Test path with one possibility
        res = graphs.findShortestMSTPath(self.mst3, "0_2", ["0_0"])
        expected = ['0_2', '0_1', '0_0']
        self.assertEqual(expected, res)

        # Test path with one possibility
        res = graphs.findShortestMSTPath(self.mst3, "0_2", ["0_1"])
        expected = ['0_2', '0_1']
        self.assertEqual(expected, res)

        # Test path with one possibility
        res = graphs.findShortestMSTPath(self.mst3, "0_2", ["0_5"])
        expected = ['0_2', '0_1', '0_4', '0_5']
        self.assertEqual(expected, res)

        # Test path with two possibilities, the shorter one is 0_0 and should
        # be chosen
        res = graphs.findShortestMSTPath(self.mst3, "0_2", ["0_5", "0_0"])
        expected = ['0_2', '0_1', '0_0']
        self.assertEqual(expected, res)

if __name__ == '__main__':
    unittest.main()

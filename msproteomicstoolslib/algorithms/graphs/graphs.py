#!/usr/bin/env python
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

def getAdjacencyList(tree):
    """
    Convert a tree into a adjacency list

    Args:
        tree(list(tuple)): a tree represented as list of edges (for example [('0', '1'), ('1', '2')] ).
    """
    adj_list = {}
    for e1, e2 in tree:
        e1l = adj_list.get(e1, [])
        e1l.append(e2)
        adj_list[e1] = e1l
        e2l = adj_list.get(e2, [])
        e2l.append(e1)
        adj_list[e2] = e2l
    return adj_list

def doBFS(tree, start_node):
    """
    Perform breadth-first-search (BFS) on a given tree

    Args:
        tree(list(tuple)): a tree represented as list of edges (for example [('0', '1'), ('1', '2')] ).
        start_node(str): starting node

    Yields:
        node(str): current node during search
    """

    adj_list = getAdjacencyList(tree)

    # Do BFS
    current_front = [ start_node ]
    already_visited = set([])
    next_front = set( current_front )
    while len(next_front) > 0:
        # print " === current front", current_front
        for node in sorted(current_front):
            # print "look at", node, "connecting: ", adj_list[node]
            next_front.update( adj_list[node] )
            yield node

        # print "next front", next_front
        already_visited.update( set(current_front) )
        next_front = set(next_front) - already_visited
        # print "next front", next_front
        current_front = set(next_front)

def findOnePath(graph, start, end, path=[]):
    """
    Finds a path in a graph between start and end. 
    """
    path = path + [start]
    if start == end:
        return path
    if not start in graph:
        return None
    for node in graph[start]:
        if node not in path:
            newpath = findOnePath(graph, node, end, path)
            if newpath: return newpath
    return None

def findShortestMSTPath(graph, start, end):
    """
    Finds a path in an MST from start to one of the end elements

    The algorithm will look for the shortest path in a minimum spanning tree
    (MST) to one of the elements contained in end. It will do a breadth-first
    search (BFS) through the MST to find the first element in "end" which has
    minimal distance to start. If there are multiple elements in "end" with
    equal distance, whichever comes first in the BFS will be chosen.

    It will then return the path between this element and the start element. 

    Args:
        graph(list(tuple)): a graph represented as list of edges (for example [('0', '1'), ('1', '2')] ).
        start(str): Starting node
        end(list(str)): List of possible end nodes (closest will be chosen)

    Returns:
        path (list(str)) : A path represented as list of nodes to be visited in that order

    """
    assert isinstance(end, list)

    for node in doBFS(graph, start):
        if node in end:
            break

    # Since the graph is an MST, any path will be good enough (there is only
    # one path between two nodes)
    return findOnePath( getAdjacencyList(graph), start, node)



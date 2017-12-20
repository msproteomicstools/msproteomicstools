# distutils: language = c++
# cython: c_string_encoding=ascii  # for cython>=0.19
# encoding: latin-1
cimport cython
cimport libc.stdlib
cimport numpy as np

cdef cppclass c_node:
    # """
    # A node in the C++ MST tree structure

    # Each node knows about its neighbors and its own id
    # """

    libcpp_string internal_id
    libcpp_vector[c_node*] neighbors

cdef cppclass cpp_tree:
    # """
    # A C++ MST tree which contains references to all nodes
    # """

    libcpp_map[libcpp_string, c_node*] nodes

cdef visit_tree_simple(cpp_tree tree, c_node * current, libcpp_unordered_set[libcpp_string] visited):
    """
    Visit all nodes in the MST
    """

    # mark as visited
    visited.insert( deref(current).internal_id )
    cdef c_node* node1
    cdef c_node* node2
    cdef libcpp_string e1
    cdef libcpp_vector[c_node*].iterator n_it = deref(current).neighbors.begin()
    while (n_it != deref(current).neighbors.end()):
        if visited.find( deref(deref(n_it)).internal_id) != visited.end():
            # have visited this neighbor already
            inc(n_it)
        else:
            # have not visited this edge (current, neighbor) yet
            node1 = current
            node2 = deref(n_it)
            # do work
            print ("visit!", deref(node1).internal_id, deref(node2).internal_id)
            # recursive call
            visit_tree_simple(tree, node2, visited)
            inc(n_it)


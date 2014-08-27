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


# from  http://stackoverflow.com/questions/12151182/python-precondition-postcondition-for-member-function-how
import functools
def condition(pre_condition=None, post_condition=None, class_invariant=None):

    # Dummy decorator (calls original function)
    def dummyDecorator(func):
        @functools.wraps(func) # preserves name, docstring, etc
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper

    # Real check of pre/post condition and class invariants
    def decorator(func):

        @functools.wraps(func) # preserves name, docstring, etc
        def wrapper(*args, **kwargs): #NOTE: no self

            # Check class invariant and pre-condition before
            if class_invariant is not None:
                assert class_invariant(args[0]), "Class invariant failed: %s" % class_invariant.__name__
            if pre_condition is not None:
                assert pre_condition(*args, **kwargs), "Pre condition failed: %s" % pre_condition.__name__

            retval = func(*args, **kwargs) # call original function or method

            # Check class invariant and post-condition after
            if post_condition is not None:
                assert post_condition(retval), "Post condition failed: %s" % post_condition.__name__
            if class_invariant is not None:
                assert class_invariant(args[0]), "Class invariant failed: %s" % class_invariant.__name__

            # Return retval, end of wrapper
            return retval
        return wrapper

    # Enable/disable checking of pre/post conditions depending on whether
    # Python assertions are turned on or off (there is no point to check all
    # conditions if assertions are not on). Returning the dummy decorator
    # should also remove some of the overhead.
    try:
        assert(False)
        # Assertions are off! Lets use the dummy decorator
        return dummyDecorator
    except AssertionError:
        return decorator

def pre_condition(check):
    return condition(pre_condition=check)

def post_condition(check):
    return condition(post_condition=check)

def class_invariant(check):
    return condition(class_invariant=check)

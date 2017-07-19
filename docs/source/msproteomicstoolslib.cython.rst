Optimized Cython 
================

This document contains information about the optimized Cython data structures
and algorithms available within this package.

The main data structures are  (see :any:`msproteomicstoolslib.data_structures` for a detailed description):

- :class:`.CyPeakgroupWrapperOnly` (wraps a C++ peakgroup)
- :class:`.CyPrecursorWrapperOnly` (wraps a C++ precursor)
- :class:`.CyPrecursorGroup` (a precursor group, a Cython version of :class:`.PrecursorGroup`)

For linear interpolation and retention time transformation, the
:class:`.CyLinearInterpolateWrapper` wraps a C++ interpolation function,
allowing access to a very fast linear interpolator. The
:class:`.CyLightTransformationData` is a Cython version of
:class:`.LightTransformationData`.

The :class:`.CyDataCacher` is a Cython class that holds pairwise RT alignment
data cached for later use.

Finally, the MST alignment algorithm can be called through
:any:`static_cy_alignBestCluster`.


Peakgroup (optimized)
^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: msproteomicstoolslib.cython._optimized.CyPeakgroupWrapperOnly
    :members:
    :undoc-members:
    :show-inheritance:

Precursor (optimized)
^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: msproteomicstoolslib.cython._optimized.CyPrecursorWrapperOnly
    :members:
    :undoc-members:
    :show-inheritance:

PrecursorGroup (optimized)
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: msproteomicstoolslib.cython._optimized.CyPrecursorGroup
    :members:
    :undoc-members:
    :show-inheritance:

CyLinearInterpolateWrapper (optimized)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: msproteomicstoolslib.cython._optimized.CyLinearInterpolateWrapper
    :members:
    :undoc-members:
    :show-inheritance:

CyLightTransformationData (optimized)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: msproteomicstoolslib.cython._optimized.CyLightTransformationData
    :members:
    :undoc-members:
    :show-inheritance:


DataCacher (optimized)
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: msproteomicstoolslib.algorithms.alignment.DataCacher.CyDataCacher
    :members:
    :undoc-members:
    :show-inheritance:



Algorithms - Alignment
======================

:mod:`MST Alignment` Module
-------------------

MST Alignment
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.algorithms.alignment.AlignmentMST.TreeConsensusAlignment
    :members:
    :undoc-members:
    :show-inheritance:
    :special-members:

.. autofunction:: msproteomicstoolslib.algorithms.alignment.AlignmentMST.getDistanceMatrix

Simple Alignment Algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm.Cluster
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.algorithms.alignment.AlignmentAlgorithm.AlignmentAlgorithm
    :members:
    :undoc-members:
    :show-inheritance:

SplineAligner
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.algorithms.alignment.SplineAligner.SplineAligner
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.algorithms.alignment.SplineAligner.TransformationError
    :members:
    :undoc-members:
    :show-inheritance:

BorderIntegration 
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: msproteomicstoolslib.algorithms.alignment.BorderIntegration.integrationBorderShortestPath
.. autofunction:: msproteomicstoolslib.algorithms.alignment.BorderIntegration.integrationBorderShortestDistance
.. autofunction:: msproteomicstoolslib.algorithms.alignment.BorderIntegration.integrationBorderReference

:mod:`Alignment Datastructures` Module
-------------------

Multipeptide 
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.algorithms.alignment.Multipeptide.Multipeptide
    :members:
    :undoc-members:
    :show-inheritance:

MRExperiment 
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.algorithms.alignment.MRExperiment.MRExperiment
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`Alignment Helper` Module
-------------------

FDREstimation
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.algorithms.alignment.FDRParameterEstimation.ParamEst
    :members:
    :undoc-members:
    :show-inheritance:

.. autofunction:: msproteomicstoolslib.algorithms.alignment.AlignmentHelper.write_out_matrix_file
.. autofunction:: msproteomicstoolslib.algorithms.alignment.AlignmentHelper.addDataToTrafo

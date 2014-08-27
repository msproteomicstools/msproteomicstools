Alignment executables 
======================

:mod:`FeatureAlignment` Module
-------------------

.. autoclass:: feature_alignment.AlignmentStatistics
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass:: feature_alignment.Experiment
    :members:
    :undoc-members:
    :show-inheritance:

.. autofunction:: feature_alignment.estimate_aligned_fdr_cutoff
.. autofunction:: feature_alignment.doMSTAlignment
.. autofunction:: feature_alignment.doParameterEstimation
.. autofunction:: feature_alignment.doReferenceAlignment
.. autofunction:: feature_alignment.main


:mod:`Noise imputation` Module
-------------------


Analysis functions
^^^^^^^^^^^^^^^
.. autofunction:: requantAlignedValues.runSingleFileImputation
.. autofunction:: requantAlignedValues.runImputeValues
.. autofunction:: requantAlignedValues.analyze_multipeptides
.. autofunction:: requantAlignedValues.analyze_multipeptide_cluster
.. autofunction:: requantAlignedValues.integrate_chromatogram
.. autofunction:: requantAlignedValues.write_out
.. autofunction:: requantAlignedValues.main


Data Structures
^^^^^^^^^^^^^^^

.. autoclass:: requantAlignedValues.ImputeValuesHelper
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: requantAlignedValues.SwathChromatogramRun
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: requantAlignedValues.SwathChromatogramCollection
    :members:
    :undoc-members:
    :show-inheritance:

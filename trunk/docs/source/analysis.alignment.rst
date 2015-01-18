Alignment executables 
======================

:mod:`FeatureAlignment` executable
-------------------

The Feature Alignment executable can be run as ::

  python feature_alignment.py 

and the for help please use ::

  python feature_alignment.py --help 

Some of the most used options are the following

fdr_cutoff
""""""""""
This is the seeding score cutoff, if a precursor has an identification in one
run with at least this score, it will be included for alignment.

max_fdr_quality
""""""""""""""""
This is the extension score cutoff. During each step of the algorithm, a
peakgroup from a new run is added to the initial seed (see above). Only if the
additional peakgroup in the new run has a score better than max_fdr_quality
will it be included in the final result.

target_fdr
""""""""""
Experimental option for dynamic parameter estimation of the fdr_cutoff
parameter. If you want to use this, please turn off fdr_cutoff (but
max_fdr_quality still needs to be set).

method
""""""
Defines the method to use for the clustering. Available options are 

* best_overall
* best_cluster_score 
* global_best_cluster_score
* global_best_overall
* LocalMST
*  LocalMSTAllCluster

Note that the MST options will perform a local, MST guided alignment while the
other options will use a reference-guided alignment. The global option will
also move peaks which are below the selected FDR threshold (while the
best_overall and best_cluster_score will not touch any peak that is below fdr_cutoff).

realign_method
""""""""""""""
Method to use to re-align retention times between pairs of runs. The following options are available:


* None: use the raw RT from the file (not recommended) 
* diRT: use only deltaiRT from the input file
* linear: perform a linear regression using best peakgroups
* splineR: perform a spline fit using R (this feature relies on the rpy2 package)
* splineR_external: perform a spline fit using R (start an R process using the command line, not tested under Windows)
* splinePy: use Python native spline from scikits.datasmooth (not recommended, very slow)
* nonCVSpline, CVSpline: splines with and without cross-validation from scipy.interpolate
* lowess: use Robust locally weighted regression (lowess smoother)
* earth : use Multivariate Adaptive Regression Splines using py-earth
* WeightedNearestNeighbour: the weighted RT of the nearest neighbours is used
* SmoothLLDMedian: a local kernel of linear differences is computed 

Recommended options are CVSpline and splineR and splineR (if you have R). Both
WeightedNearestNeighbour and SmoothLLDMedian gave acceptable results.




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

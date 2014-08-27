Format
======

Several classes and functions to deal with common mass spectrometric format (mostly dealing with File I/O).

:mod:`Transformation Collection` Module
-------------------

TransformationCollection
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.format.TransformationCollection.TransformationCollection
    :members:
    :undoc-members:
    :show-inheritance:

LightTransformationData
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.format.TransformationCollection.LightTransformationData
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`File Reader` Module
-------------------

SWATHScoringReader
^^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.format.SWATHScoringReader.ReadFilter
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.format.SWATHScoringReader.SWATHScoringReader
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.format.SWATHScoringReader.OpenSWATH_SWATHScoringReader
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass:: msproteomicstoolslib.format.SWATHScoringReader.mProphet_SWATHScoringReader
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass:: msproteomicstoolslib.format.SWATHScoringReader.Peakview_SWATHScoringReader
    :members:
    :undoc-members:
    :show-inheritance:

.. autofunction:: msproteomicstoolslib.format.SWATHScoringReader.inferMapping

:mod:`Data Matrix` Module
-------------------

Functions for handling the output data matrix

MatrixWriters
^^^^^^^^^^^^^^

.. autofunction:: msproteomicstoolslib.format.MatrixWriters.getwriter

.. autoclass:: msproteomicstoolslib.format.MatrixWriters.IWriter
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.format.MatrixWriters.CsvWriter
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.format.MatrixWriters.XlsWriter
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.format.MatrixWriters.XlsxWriter
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`Spectral library` Module
-------------------

Functions for handling SpectraST spectral library format

Spectral library handler
^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.format.speclib_db_lib.Library
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.format.speclib_db_lib.SequenceHandler
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: msproteomicstoolslib.format.speclib_db_lib.Spectra
    :members:
    :undoc-members:
    :show-inheritance:

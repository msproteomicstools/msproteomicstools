DataStructures - Alignment
==========================

This document contains information about the data structures used in the TRIC algorithm.

- :class:`.Run` contains all data pertaining to a LC-MS/MS run, particularly references to measured precursors
- :class:`.PrecursorGroup` represents a set of precursors (e.g. precursors deriving from the same peptide sequence but identified by different charge states and isotopic labelling); see also :class:`.CyPrecursorGroup` for a Cython implementation
- A Precursor represents a single precursor (e.g. a single measured analyte with a precursor m/z identified by its chemical formula, charge state and isotopic labelling)
    - :class:`.PrecursorBase` is a base implementation of a Precursor
    - :class:`.Precursor` is the implementation of a Precursor using minimal memory, see also :class:`.CyPrecursorWrapperOnly` for a Cython implementation
    - :class:`.GeneralPrecursor` is the default implementation of a Precursor
- A peak group represents a single RT region in the chromatogram of a single Precursor
    - :class:`.PeakGroupBase` is a base implementation of a Peakgroup
    - :class:`.MinimalPeakGroup` is the implementation of a Peakgroup using minimal memory, see also :class:`.CyPeakgroupWrapperOnly` for a Cython implementation
    - :class:`.GeneralPeakGroup` is the default implementation of a Peakgroup
    - :class:`.GuiPeakGroup` is the implementation used by the GUI


:mod:`Run` Module
-------------------

Run
^^^

.. autoclass:: msproteomicstoolslib.data_structures.Run.Run
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`PrecursorGroup` Module
----------------------------

PrecursorGroup
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.PrecursorGroup.PrecursorGroup
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`Precursor` Module
------------------------

PrecursorBase
^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.Precursor.PrecursorBase
    :members:
    :undoc-members:
    :show-inheritance:

GeneralPrecursor
^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.Precursor.GeneralPrecursor
    :members:
    :undoc-members:
    :show-inheritance:

Precursor
^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.Precursor.Precursor
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`PeakGroup` Module
-----------------------

PeakGroupBase
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.PeakGroup.PeakGroupBase
    :members:
    :undoc-members:
    :show-inheritance:

MinimalPeakGroup
^^^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.PeakGroup.MinimalPeakGroup
    :members:
    :undoc-members:
    :show-inheritance:

GuiPeakGroup
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.PeakGroup.GuiPeakGroup
    :members:
    :undoc-members:
    :show-inheritance:

GeneralPeakGroup
^^^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.PeakGroup.GeneralPeakGroup
    :members:
    :undoc-members:
    :show-inheritance:

DataStructures - Basic
=======================

:mod:`Aminoacides` Module
----------------------------

Aminoacid
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.aminoacides.Aminoacid
    :members:
    :undoc-members:
    :show-inheritance:


Aminoacides
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.aminoacides.Aminoacides
    :members:
    :undoc-members:
    :show-inheritance:


:mod:`Modifications` Module
----------------------------

Modification
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.modifications.Modification
    :members:
    :undoc-members:
    :show-inheritance:


Modifications
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.modifications.Modifications
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`Peak` Module
-------------------

Peak
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.peak.Peak
    :members:
    :undoc-members:
    :show-inheritance:


:mod:`Peptide` Module
----------------------

Peptide
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.peptide.Peptide
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`Residues` Module
----------------------

Residues
^^^^^^^^^^^^^^

.. autoclass:: msproteomicstoolslib.data_structures.Residues.Residues
    :show-inheritance:

:mod:`DDB` Module
-------------------

DDB
^^^^^^^^^^^^^^

Abstraction layer to the `2DDB <http://sourceforge.net/projects/twoddb/>`_  software framework.


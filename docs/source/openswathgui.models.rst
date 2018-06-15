OpenSWATH - GUI Models
============

The main models used by the GUI are the :class:`.PeptideTree` and the
:class:`MSData.DataModel <.DataModel>` models. Internally, :class:`.PeptideTree`
uses :class:`.ChromatogramTransition` to store access to single rows in the
tree data structure while :class:`MSData.DataModel <.DataModel>` uses
:class:`.SwathRunCollection` to keep track of multiple SWATH-MS runs.

:mod:`MSData Data model` Module
-------------------

Contains classes that provide access to the raw data

MSData
^^^^^^

.. autoclass:: openswathgui.models.MSData.DataModel
    :members:
    :undoc-members:
    :show-inheritance:

.. 
    PrecursorModel
    ^^^^^^^^^^^^^^

    .. autoclass:: openswathgui.models.MSData.PrecursorModel
        :members:
        :undoc-members:
        :show-inheritance:


:mod:`TreeModels` Module
-------------------

Contains classes that provide access to the hierarhical tree container protein,
precursor, peptide and transition level data.

While :class:`.TreeNode` and :class:`.TreeModel` are generic models for trees
and nodes, the derived classes :class:`.PeptideTreeNode` and
:class:`.PeptideTree` are implementations specific to TAPIR.

TreeNode
^^^^^^^^^^^^^^

.. autoclass:: openswathgui.models.TreeModels.TreeNode
    :undoc-members:
    :show-inheritance:

TreeModel
^^^^^^^^^^^^^^

.. autoclass:: openswathgui.models.TreeModels.TreeModel
    :undoc-members:
    :show-inheritance:

PeptideTreeNode
^^^^^^^^^^^^^^

.. autoclass:: openswathgui.models.PeptideTree.PeptideTreeNode
    :members:
    :show-inheritance:

PeptideTree
^^^^^^^^^^^^^^

.. autoclass:: openswathgui.models.PeptideTree.PeptideTree
    :members:
    :undoc-members:
    :show-inheritance:


:mod:`SWATH MS Run` Module
-------------------

Raw chromatographic data is handled using the :class:`.SwathRunCollection` class which can either hold references to mzML or to SqMass data.

SwathRunCollection
^^^^^^^^^^^^^^^^^^

.. autoclass:: openswathgui.models.SwathRunCollection.SwathRunCollection
    :members:
    :undoc-members:
    :show-inheritance:

SqMass
^^^^^^

.. autoclass:: openswathgui.models.SqlSwathRun.SqlSwathRun
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: openswathgui.models.SqlDataAccess.SqlDataAccess
    :members:
    :undoc-members:
    :show-inheritance:

MzML File
^^^^^^^^^

.. autoclass:: openswathgui.models.SwathRun.SwathRun
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: openswathgui.models.SingleChromatogramFile.SingleChromatogramFile
    :members:
    :undoc-members:
    :show-inheritance:



:mod:`ChromatogramTransition` Module
-------------------

ChromatogramTransition
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: openswathgui.models.ChromatogramTransition.ChromatogramTransition
    :members:
    :undoc-members:
    :show-inheritance:


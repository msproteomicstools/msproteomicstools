OpenSWATH - GUI Views
============

:mod:`PeptideTreeView` Module
-------------------

PeptideTreeView
^^^^^^^^^^^^^^

.. autoclass:: openswathgui.views.PeptideTree.PeptidesTreeView
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`Plot` Module
-------------------

Plot
^^^^^^^^^^^^^^

There are two implementations of the plotting view, one uses guiqwt (may not be
present on all systems) and the other one uses plain Qwt (should be safer).


.. autoclass:: openswathgui.views.Plot.CurveItemModel
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: openswathgui.views.Plot.GuiQwtMultiLinePlot
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: openswathgui.views.Plot.QwtMultiLinePlot
    :members:
    :undoc-members:
    :show-inheritance:

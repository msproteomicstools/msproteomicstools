TAPIR GUI executable
============================

This documents the functions in the TAPIR.py executable and also the underlying
views and models (see TOC at the end). The main window is :class:`.MainWindow`
which contains a single main widget :class:`.ApplicationView`. The
ApplicationView is the main widget, it contains the left side tree structure
(peptide tree) and the right side graph area.

TAPIR Main Window 
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: TAPIR.MainWindow
    :members:
    :undoc-members:
    :show-inheritance:
    :exclude-members: center,initUI


TAPIR ApplicationView (central window) 
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: TAPIR.ApplicationView
    :undoc-members:
    :show-inheritance:
    :members:

TAPIR Widgets
^^^^^^^^^^^^^^

.. autoclass:: TAPIR.GraphArea
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: TAPIR.PeptideTreeWidget
    :members:
    :undoc-members:
    :show-inheritance:


TAPIR Config Dialog
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: TAPIR.ConfigDialog
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: TAPIR.Settings
    :members:
    :undoc-members:
    :show-inheritance:

TAPIR Models and Views
============================

.. toctree::

    openswathgui.models
    openswathgui.views


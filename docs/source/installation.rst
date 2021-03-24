.. include:: ./substitutions.rst

Installing brutifus
===================

brutifus is available on pypi, which makes its installation easier than ever.
In a terminal, type:
::

   pip install brutifus

And that should take care of things.

The most recent release of brutifus is also available for download from its
`Github repository <https://github.com/fpavogt/fcmaker/releases/latest/>`_.
Interested users can fork the ``develop`` branch of the brutifus repository if they want to get access
to the latest updates not yet released. *Push requests* for bug fixes and new features are
welcome and will be examined in detail.

Requirements
------------
brutifus is written in Python 3. It is compatible with he following versions:

.. literalinclude:: ../../setup.py
    :language: python
    :lines: 19

brutifus also relies on the following python packages, that will get automatically installed by pip:

.. literalinclude:: ../../setup.py
    :language: python
    :lines: 20-30

Testing the installation
------------------------

In a terminal shell, try calling the high-level brutifus help:

.. code-block:: none

   brutifus --help

This should return the following information:

   .. literalinclude:: brutifus_help_msg.txt
       :language: none

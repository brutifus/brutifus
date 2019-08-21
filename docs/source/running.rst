
Running brutifus
================

The spirit of brutifus is that each user can choose, depending on the object at hand & the
quality of the data, what processing steps are warranted. These are governed by the
``procsteps_brutifus.yaml`` file. In parallel, the user can define all relevant parameters using
the ``params_brutifus.yaml`` file. The two files (that rely on YAML syntax) are `shipped as 
supplementary files with the code <https://github.com/fpavogt/brutifus/blob/master/brutifus/exec_scripts/>`_. 
The following high-level entry point will create local copies in your current location (assumed to be your
favorite processing area). In a terminal: :: 

    cd ~/where/ever/you/want
    brutifus --setup


Doing so will create a local copy of both files in the current directory, together with generic
working directories.

.. note:: 

    I very much suggest to run a ``brutifus --setup`` for each project, rather than using 
    a single processing area for all projects.


Setting parameters with ``params_brutifus.yaml``
------------------------------------------------

All the brutifus parameters with a `scientific impact` can be specified inside this file.

.. literalinclude:: ../../brutifus/exec_scripts/params_brutifus.yaml
   :language: yaml
   :linenos:

.. todo::
   
    A detailed list of what everything is/does would be nice. Until then, look at the 
    comments inside the file for a clarification of what everything does.


Setting the processing steps with ``procsteps_brutifus.yaml``
-------------------------------------------------------------

brutifus is designed to batch process all the individual spectra inside a given IFU datacube 
individually, exploiting multiple cpus (when available) to gain speed. 

The different steps to be executed, in sequential order, are defined in ``procsteps_brutifus.py``. 
Each step can be re-arranged as needed by the user. In practice, it is mostly the plotting steps
that will be most duplicated to visualize the different intermediate products.

Each step contains a handful of common parameters, including a step name, whether to run it (True) or
not (False), and the suffix added to the associated product file. Arguments (such as cut levels for 
images) are passed via the ``args`` container. Of note are the ``name_in`` and ``name_out`` values:
these allow the user to give `short reference tags` to specific intermediate products, in order to
easily feed them to subsequent steps.



.. literalinclude:: ../../brutifus/exec_scripts/procsteps_brutifus.yaml
   :language: yaml
   :linenos:

.. todo::
   
    A detailed list of what everything is/does would be nice. Until then, look at the 
    comments inside the file for a clarification of what everything does.



Launching the post-processing 
-----------------------------

Having prepared both ``params_brutifus.yaml`` and ``procsteps_brutifus.yaml``, you are ready to launch
brutifus. This is done via the following high-level entry point::

    >>> brutifus -e procsteps_brutifus.yaml params_brutifus.yaml
    
If all goes according to plan, brutifus will start processing your data, one step after the 
other, printing some info along the way. **Some steps can take time (up to several hours) !** 
Just be patient - progress indicators will help you monitor each task, if you set ``'verbose': True`` 
inside ``params_brutifus.yaml``. 





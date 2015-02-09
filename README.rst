=======================
Loop modeling benchmark
=======================

(...Description...)

This benchmark is part of the `Modeling protocol benchmarks and protocol captures <https://kortemmelab.ucsf.edu/benchmarks>`_ website.

---------
Licensing
---------

The contents of the repository *where possible* are licensed under the Creative Commons Attribution 4.0 International License. The license only applies to files which either: i) include the license statement; or ii) which are explicitly listed in some file in the repository as being covered by the license. All other files may be covered under a separate license. The LICENSE file in the root of this repository is present only for the convenience of the user to indicate the license which covers any novel content presented herein.

This protocol capture is based off the original captures from the Smith & Kortemme papers listed above however most of
the output directories have been excluded here to reduce the size of the repository.

The original output directories can be found in the `RosettaCommons repositories <https://github.com/RosettaCommons/demos/tree/master/protocol_capture/2010/backrub_seqtol>`_ or at http://kortemmelab.ucsf.edu/data/.

-------------------------
Downloading the benchmark
-------------------------

The benchmark is hosted on GitHub. The most recent version can be checked out using the `git <http://git-scm.com/>`_ command-line tool:

::

  git clone https://github.com/Kortemme-Lab/loop_modeling.git

---------------------------
Directories in this archive
---------------------------

This archive contains the following directories:

- *data* : contains some data from the Smith and Kortemme 2010 publication about the input structures and proteins.
- *input* : contains the input files for the benchmark.
- *output* : these directories are empty by default. This is the default output location for protocols if they are run on the local machine.
- *output/sample* : contains sample output data that can be used to test the analysis script.
- *analysis* : contains the analysis scripts used to analyze the output of a prediction run. All protocols are expected to produce output that will work with the analysis scripts.
- *protocols* : contains the scripts needed to run a job. The scripts for a protocol are provided in a specific subdirectory.
- *hpc* : contains scripts that can be used to run the entire benchmark using specific cluster architectures. For practical reasons, a limited number of cluster systems are supported. Please feel free to provide scripts which run the benchmark for your particular cluster system.

--------------------------------------
Analysis
--------------------------------------

The analysis script generates ...some_number... metrics which can be used to evaluate the results of the loop modeling simulations.
The analysis scripts are described in more detail in analysis/README.rst.

===================================================
Protocol 1: Next-Generation Kinematic Closure (NGK)
===================================================

Created by: ...authors...

Software suite: Rosetta

Protocol directory: rosetta

============================================
Protocol 2: Kinematic Closure (KIC) - Legacy
============================================

-------------------
General Information
-------------------

Created by: ...authors...

Software suite: Rosetta

Protocol directory: rosetta

-------------------
Description
-------------------

...todo...

-------------------
Instructions
-------------------

...todo...

------------------------
Protocol capture scripts
------------------------

...todo...

-------------------
Common Flags
-------------------

_____________
General flags
_____________

+----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------+
+============================+===========================================================================================================================================+
| -s 	                     | This flag specifies the starting structure.                                                                                               |
+----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------+



___________________
Loop modeling flags
___________________



+---------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
+===========================+===================================================================================================================================================================+
+---------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ...some_option...         | ...some_description...                                                                            |
+---------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+


----------------------
Example command lines
----------------------

::

  rosetta/source/bin/loopmodel.linuxgccrelease -database rosetta/database
  -s ...todo...

----------------------------
Supporting tool versions
----------------------------

This protocol capture has been tested with:

- Python 2.6.6
- Python 2.7.8

-------------------------------------------------
References to published works using this protocol
-------------------------------------------------

...todo...

Try to follow the existing format e.g.
Smith, CA, Kortemme, T. Structure-Based Prediction of the Peptide Sequence Space Recognized by Natural and Synthetic PDZ Domains. 2010. J Mol Biol 402(2):460-74. `doi: 10.1016/j.jmb.2010.07.032 <http://dx.doi.org/10.1016/j.jmb.2010.07.032>`_.




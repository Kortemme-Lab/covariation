====================================
Sequence profile and covariation benchmark
====================================

Evolutionary pressures on protein structure and function have shaped the amino acid sequences of today's naturally occurring proteins. Consequently, the sequences of natural proteins are nearly optimal for their structures. Natural protein sequences therefore provide valuable information for evaluating the accuracy of computational protein design. The purpose of this benchmark is to evaluate the extent to which protein design can recapitulate properties of naturally occurring proteins, including amino acid sequence preferences (“sequence profiles”) and patterns of amino acid covariation.

This benchmark includes:

- todo:

---------
Licensing
---------

The contents of the repository *where possible* are licensed under the MIT License. The license only applies to files which either: i) include the license statement; or ii) which are explicitly listed in some file in the repository as being covered by the license. All other files may be covered under a separate license. The LICENSE file in the root of this repository is present only for the convenience of the user to indicate the license which covers any novel content presented herein.

-------------------------
Downloading the benchmark
-------------------------

The benchmark is hosted on GitHub. The most recent version can be checked out using the `git <http://git-scm.com/>`_ command-line tool:

::

  git clone https://github.com/Kortemme-Lab/covariation.git

---------------------------
Directories in this archive
---------------------------

This archive contains the following directories:

- *input* : contains the input files for the benchmark. Input files specific to a particular protocol are in a subdirectory named after the protocol. The input files are described in more detail in input/README.rst.
- *output* : these directories are empty by default. This is the default output location for protocols if they are run on the local machine.
- *output/sample* : contains sample output data that can be used to test the analysis script.
- *analysis* : contains the analysis scripts used to analyze the output of a prediction run. All protocols are expected to produce output that will work with the analysis scripts.
- *protocols* : contains the scripts needed to run a job. The scripts for a protocol are provided in a specific subdirectory.
- *hpc* : contains scripts that can be used to run the entire benchmark using specific cluster architectures. For practical reasons, a limited number of cluster systems are supported. Please feel free to provide scripts which run the benchmark for your particular cluster system.

=========
Protocols
=========

This repository contains one protocol which can be used to run the benchmark. We welcome the inclusion of more protocols.
Please contact support@kortemmelab.ucsf if you wish to contribute towards the repository.

Each protocol is accompanied by specific documentation in its protocol directory.

--------------------------------------
Protocol 1: Protocol name? todo
--------------------------------------

Created by: Noah Ollikainen

Software suite: Rosetta

Protocol directory: protocols/todo

__________
References
__________

Computational protein design quantifies structural constraints on amino acid covariation. 2013.
Ollikainen N, Kortemme T. PLoS Comput Biol 9(11):e1003313. `doi: 10.1371/journal.pcbi.1003313 <http://dx.doi.org/10.1371/journal.pcbi.1003313>`_. Epub 2013 Nov 14.

========
Analysis
========

The same set of analysis scripts is used by all protocols. Conceptually, the analysis scripts should be a black box that
is separated from the output of each protocol by an interface.

The analysis scripts generates ??? metrics which can be used to evaluate the results of the covariation simulations (as
well as a series of plots???). The scripts are described in more detail in analysis/README.rst.

=============================
Generic cluster command lines
=============================

The command lines in this directory demonstrate how to run the Rosetta protocols (protocols/rosetta) on a computational cluster or high-performance computing resource.

The main distinction between these scripts and those in the protocols folder is that the -nstruct parameter (number of generated structures) is set to 1 rather than 500 i.e. it is intended that the scripts here should be run 500 times in parallel whereas the protocols/rosetta scripts are intended for a single processor system which creates 500 structures in sequence. Running the Rosetta simulations in parallel involves some bookkeeping with respect to the uniqueness of filenames. To this end, the Rosetta -out::suffix parameter is used to ensure that the output of parallel jobs do not overwrite each other. The scripts here use Sun Grid Engine environment variables ($SGE_TASK_ID) to demonstrate where the appropriate changes should be made for your particular cluster system.


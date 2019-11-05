
================
SigProfilerJulia
================

|

SigProfilerJulia aims to extract the mutational processes active in a tumor cohort. It is a fast implementation in Julia of the original SigProfiler written in MATLAB.
The original motivation for this implementation was that we did not have the MATLAB license to run SigProfiler with more than 12 cores in our local cluster.

This software was developed by Oriol Pich, Ferran Muiños and Jordi Deu-Pons and was used in `Pich et al, The mutational footprints of cancer therapies <https://www.biorxiv.org/content/10.1101/683268v1>`_.

|

|

----------------------
Using SigProfilerJulia
----------------------

The input to SigProfilerJulia consists of a tab-separated values (tsv) file containing the counts of each mutational type or channel each of the tumor samples.
The amount and type of channels depend on the nature of the variants (e.g. single-nucleotide variants, double-substitutions or short-indels). A single-nucleotide variant extraction typically relies on 96 channels (the trinucleotide context of the variant plus the alternate allele).

We currently provide one script to create the 96-channel single-nucleotide variant count matrix in a easy and fast way. Please see further details and one example in the `preprocessing folder  <https://bitbucket.org/bbglab/sigprofilerjulia/src/master/preprocessing/>`_.
Double-substitutions or indel channel creation is not yet implemented here, however please refer to `The mutational footprints of cancer therapies repository  <https://bitbucket.org/bbglab/mutfootprints/src/master/>`_ where the assembly of the corresponding counts matrices is implemented (and also how to postprocess the outputs).

We also provide two postprocessing scripts, which will help the user to decide the number of signatures active in the tumors, as well
as finding similaries with signatures in COSMIC and plotting all the extracted profiles. Please see further details in the `postprocessing folder  <https://bitbucket.org/bbglab/sigprofilerjulia/src/master/postprocessing/>`_.

|

---------------------------------------------
Software dependencies to run SigProfilerJulia
---------------------------------------------

Julia installables must be downloaded from the `Julia downloads site  <https://julialang.org/downloads/>`_.

This package has been developed and tested using Julia 1.0.0.

Several Julia packages are needed for SigProfilerJulia to run properly. They can be found inside julia_dependencies.jl.
To install them, please run:

.. code-block::

  $ julia julia_dependencies.jl

Where ``julia`` is the Julia binary. The path will vary depending on the OS system. Please see the `Platform specific instructions <https://julialang.org/downloads/platform.html>`_ for more insight.

|

---------------------------
How to run SigProfilerJulia
---------------------------

SigProfilerJulia is run in the command line. For a specific matrix count, it will explore a range of numbers of active mutational processes in the samples.
SigProfilerJulia also has a random seed implemented, to guarantee reproducibility between runs.

|

Example
-------

The following command carries out mutational signatures extraction in 21 WGS tumor breast samples.
We will explore the presence from 2 to 10 mutational signatures, using 100 iterations (we strongly recommend using more than 1000 in real cases) and 3 CPUs.

.. code-block::

  $ julia SigProfilerJulia.jl -f preprocessing/21_breast_WGS.snvs.txt -s 2 -m 10 -i 100 --workers 3 --outpath test/output

Where ``julia`` is the Julia executable.

|

Arguments
---------

* ``--mutation_file`` [alias= ``-f`` ]

  - File with the mutation count, where the sample names are in the header.

* ``--min_signatures_to_extract`` [alias= ``-s`` ]

  - The minimum number of signatures to extract.

* ``--max_signatures_to_extract`` [alias= ``-m`` ]

  - The maximum number of signatures to extract.

* ``--outpath`` [alias= ``-o`` ]

  - Path where the output files will be stored.

* ``--iterations`` [alias= ``-i``, default=1024 ]

  - Total NMF iterations performed. It is strongly recommended to run more than 1000.

* ``--workers`` [alias= ``-n``, default=1 ]

  - Number of local cores or slurm processes to run the code in parallel. Please note for slurm processes the ``--slurm`` flag is needed

* ``--weak`` [alias= ``-w``, default=0 ]

  - Whether to remove weak mutation types. If unsure, use default.

* ``--slurm`` [optional]

  - If you want to run SigProfilerJulia using SLURM, write this flag.

|

Outputs
---------

SigProfilerJulia produces several output files for each of the number of signatures evaluated. If SigEval is the signature we are examining, the output files would be the following:

* ``exposures_SigEval`` : matrix of exposures of SigEval signatures.

* ``processes_SigEval`` : matrix of processes of SigEval signatures.

* ``processesStabAvg_SigEval`` : the stability of each of the SigEval processes. Stability essentially reflects how consistently each mutational process is recovered when the input count data is randomized, thereby giving a sense of how reliable the mutational process is.

* ``avgStability_SigEval`` : average stability of SigEval processes.

* ``exposures_fitting_SigEval`` : matrix of exposures after removeAllSingleSignatures module is run. We recommend to use this as a matrix of exposures.

* ``reconError_SigEval`` : per sample reconstruction error of the original matrix count using the exposure count of SigEval processes using norm as a measure.

* ``reconErrorPercentage_SigEval`` : percentage of the above.

* ``simError_SigEval`` : per sample reconstruction error compared to the original matrix count using the exposure count of SigEval processes using cosine distance as a measure.

* ``avgReconstructionError_SigEval`` : average reconstruction error of the original matrix count using the exposure count of SigEval processes.

* ``avgReconstructionErrorPercentage_SigEval`` : percentage of the above.

|

---------------------------
Understanding the output
---------------------------

As mentioned before, please see further detail in the `postprocessing folder  <https://bitbucket.org/bbglab/sigprofilerjulia/src/master/postprocessing/>`_.

|

----------------------------
How to cite SigProfilerJulia
----------------------------

SigProfilerJulia has been developed by Oriol Pich, Ferran Muiños and Jordi Deu-Pons. If you have used this software, please cite:

.. admonition:: Citation
   :class: note

   Oriol Pich, Ferran Muiños, Martijn Paul Lolkema, Neeltje Steeghs, Abel Gonzalez-Perez, Nuria Lopez-Bigas, `The mutational footprints of cancer therapies <https://www.biorxiv.org/content/10.1101/683268v1>`_


The original SigProfiler implementation is described in:

.. admonition:: Citation
   :class: note

   Ludmil B. Alexandrov, Serena Nik-Zainal, David C. Wedge, Peter J. Campbell, Michael R. Stratton, `Deciphering Signatures of Mutational Processes Operative in Human Cancer <https://doi.org/10.1016/j.celrep.2012.12.008>`_

And can be found `here <https://mathworks.com/matlabcentral/fileexchange/38724-sigprofiler>`_

This software is licensed under the `3-clause BSD License <https://opensource.org/licenses/BSD-3-Clause>`_.
=======================
Preprocessing utilities
=======================

|

SigProfilerJulia needs a input matrix which contains the count of mutations falling in specific context at the sample level. We provide here
one script that will generate the 96 SNV channels (that is, the trinucleotide context and the alternate).

|

---------------------
Software dependencies
---------------------

We recommend installing python and the script dependencies using conda:

.. code-block::

  $ conda create -n sigprofilerjulia python=3.6 tqdm click pandas matplotlib=3.1.0 scipy seaborn=0.9.0
  $ conda activate sigprofilerjulia
  $ conda install -c bbglab bgreference

Alternatively, the packages can be installed using pip:

.. code-block::

  $ pip install tqdm click pandas matplotlib==3.1.0 scipy seaborn==0.9.0 bgreference

For OSX users, all the dependencies must be installed using pip. The scripts were not tested in Windows systems.


|

---------------------
How to run the script
---------------------

The script is run in the command line. The input file should be a tsv file with at least the following columns:

"CHROM": the chromosome of the variant. Please note the script will get the context for chromosomes 1 to 22 + XY.
"POS": the position of the variant.
"REF": reference allele.
"ALT": alternate allele.
"SAMPLE": name of the sample.

Otherwise the script will give an error. Using "chrX" or "X" is equivalent.

Four possible reference genomes can be used out of the box to get the context, hg19, hg38, mm9 and mm10. This must be specified while executing the script. The first time the script is run it will automatically download the genome of reference, so it will take some minutes depending on your bandwidth.
Additional genome references can be installed using the `bgreference package  <https://bitbucket.org/bgframework/bgreference/>`_ .

The script it will internally check whether your reference alleles match the corresponding alleles in the reference genome.
Be aware that only SNVs are considered INDELS and DBS will be removed (the generation of these matrices might be added in the future, meanwhile please check `The mutational footprints of cancer therapies  <https://bitbucket.org/bbglab/mutfootprints/src/master/>`_  where we implemented these features in the formatting section.

|

Example
-------
The following command will generate the matrix using the list of mutations in the whole-genomes of 21 breast carcinomas (provided in this repository), as an input. These mutations are mapped to the hg19 assembly of the human genome.

.. code-block::

  $ python create_matrix_input.py  --input_file 21_breast_WGS.txt.gz --genome_reference hg19 --output_file 21_breast_WGS.snvs.txt

|

Arguments
---------

* ``--input_file``

  - File with the mutations. "CHROM", "POS", "REF", "ALT" and "SAMPLE" should be in the header.

* ``--genome_reference``

  - The reference genome where the mutations where mapped to. Please read carefully the instructions for the reference genomes above.

* ``--output_file``
  - File with the mutation count, where the sample names are in the header.

|


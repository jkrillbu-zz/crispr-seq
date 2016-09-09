.. CRISPR-Seq documentation master file, created by
   sphinx-quickstart on Wed Sep  7 13:12:06 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CRISPR-Seq Workflow Documentation
======================================

	Display a table of contents (left) and logo (right).


What is CRISPR-Seq?
-------------------

	A sequence analysis pipeline for CRISPR/cas9 edited DNA. Starting from multiplex single-end reads and a gRNA annotation of predicted cutsites, the pipeline produces 
		1. aligned bam files per sample
		2. characterization of all indels that overlap a predicted cutsite
		3. quantification of indel reads versus total reads for each sample/gene
		3. statistical significance of indel fractions
		4. QC of indel size detection accuracy per gene
		5. indel distribution plots per gene
		6. interactive sunburst to investigate population dynamics of indels within a sample/gene

	Current design assumes long 300 bp single-end reads with predicted gRNA binding site near the center of the read.

	Sequencing the cutsite of CRISPR/cas9 experiments is useful for validating that the system is causing significant loss of function. Experiments resulting in significant LOF fractions must have sufficient Cas9 activity, guide efficiency, and tolerance to the loss-of-function mutation. Finally, the ability to track population dynamics over time lends itself to many creative applications. 

Algorithm Description
---------------------

	Algorithm description. 

Running CRISPR-Seq
------------------

.. toctree::
   :maxdepth: 2

   firecloud
   docker

Citing
------

Help
----

	Please contact mburger@broadinstitute.org with any questions.


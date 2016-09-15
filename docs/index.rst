.. CRISPR-Seq documentation master file, created by
   sphinx-quickstart on Wed Sep  7 13:12:06 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CRISPR-Seq Workflow Documentation
======================================

	Display a table of contents (left) and logo (right).

	Sequencing the cutsite of CRISPR/cas9 experiments is important validating that the system is causing significant loss of function. Experiments resulting in significant LOF fractions must have sufficient Cas9 activity, guide efficiency, and tolerance to the loss-of-function mutation. Finally, the ability to track population dynamics over time lends itself to many creative applications. 

	CRISPR-Seq is a sequence analysis pipeline for CRISPR/cas9 edited DNA.

.. toctree::
   :maxdepth: 2

   index

Why use CRISPR-Seq?
===================

	Simple inputs (multiplex single-end reads and a gRNA annotation of predicted cutsites)

	Higher accuracy indel detection with single-end reads using predicted cutsites.

	Easy for experimentalists to run using the FireCloud website (Does not require use of command line)

	Cheap (Approximately x dollars per GB of FASTQ file). Billing is through Google Cloud services.

	Source code available dependency free using Docker image for computationalists. 

	Output includes:

		- aligned bam files per sample
		- characterization of all indels that overlap a predicted cutsite
		- quantification of indel reads versus total reads for each sample/gene
		- statistical significance of indel fractions
		- QC of indel size detection accuracy per gene
		- indel distribution plots per gene
		- interactive sunburst to investigate population dynamics of indels within a sample/gene


Algorithm Description
=====================

	Current design assumes long 300 bp single-end reads with predicted gRNA binding site near the center of the read.

Running CRISPR-Seq
==================

	FireCloud is the recommended method of running CRISPR-Seq for all users who don't need to modify the workflow. For users who need tweak the execution of tasks or have a different preffered computation environment, a Docker image is available with all source code. 

.. include:: firecloud.rst


Citing
======

Help
====

	Please contact mburger@broadinstitute.org with any questions.


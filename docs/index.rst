.. CRISPR-Seq documentation master file, created by
   sphinx-quickstart on Wed Sep  7 13:12:06 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CRISPR-Seq Workflow Documentation
======================================

	Sequencing the predicted cut sites of CRISPR/cas9 experiments is an effective method of validating the CRISPR/cas9 system is creating loss-of-function (LOF) mutations. Abundant LOF allele fractions indicate sufficient Cas9 activity, guide efficiency, and tolerance to LOF mutation. Sequence analysis of predicted cut sites also facilitates studying the complex population dynamics of CRISPR/cas9 edited cells under positive or negative selection. 

	The CRISPR-Seq analysis workflow inputs single-end targeted sequencing reads that span predicted CRISPR/cas9 cut sites and outputs an analysis of LOF allele fractions and detailed indel descriptions. The CRISPR-Seq algorithm is more accurate than traditional indel callers at detecting large indels (>20bp) by using the predicted cut sites based on gRNA sequence which are unique to CRISPR/cas9 experiments. Convient options for running the analysis pipeline exist for both computational and laboratory scientists.    

Why use CRISPR-Seq?
===================

	High accuracy 

		- Improved detection of large indels (>20bp) using predicted cut sites

	Run with FireCloud

		- Easy to use web interface for experimentalists 
		- Cheap (2 GB FASTQ file costs approximately $0.40 for computation and $0.33 per month for storage)
		- Billing is managed by Google Cloud services  

	Simple inputs

		- single-end reads in FASTQ form 
		- barcode annotation (multiplex only)
		- gRNA annotation 
		- negative controls 

	Comprehensive output

		- aligned bam files per sample
		- characterization of all indels that overlap a predicted cutsite
		- quantification of indel reads versus total reads for each sample/target pair
		- statistical significance of indel allele fractions
		- QC of indel size detection accuracy per target
		- indel distribution plots per target 
		- sunburst plots to investigate population dynamics

	Open source

		- Docker image with all source code
		- Option to run workflow with Snakemake


Algorithm Description
=====================

	The algorithm is fine-tuned for detecting indels in 300nt single-end reads where the predicted gRNA binding site is near the center of the read. 

	.. image:: _static/SE_design.png
	   :scale: 70%

	|

	The alignment is a two step process. First, a basic Smith-Waterman alignment identifies wild-type reads and small indels. Second, a search for reads with fragments that have high quality mappings to the reference either before or after the 50bp region around the predicted cutsite is used to identify reads with large indels. 

	.. image:: _static/Two_Step_Algorithm.png



Running CRISPR-Seq
==================

	FireCloud is the recommended method of running CRISPR-Seq for all users who don't need to modify the workflow. For users who need tweak the execution of tasks or have a different preffered computation environment, a Docker image is available with all source code. 

.. include:: firecloud.rst


Citing
======

	Tothova Z, Krill-Burger JM, Popova KD, Landers CC, Sievers QL, Yudovich D, Belizaire R, Aster JC, Morgan EA, Tsherniak A, Ebert BL. Multiplex CRISPR/Cas9-Based Genome Editing in Human Hematopoietic Stem Cells Models Clonal Hematopoiesis and Myeloid Neoplasia. Cell Stem Cell. 2017. 21(4): 547-555. PMID:28985529

Help
====

	Please contact mburger@broadinstitute.org with any questions.


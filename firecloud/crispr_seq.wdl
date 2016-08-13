
workflow crisprSeqWorkflow {
	File barcodes_list
	File barcodes_fastq
	File reads_fastq

	File ref_idxs
	File ref_fasta

	File indel_cut_interval
	File indel_cut_site
	File indel_genes
	File neg_cntrl

	call splitBarcodesTask {
		input: 
		barcodes_list = barcodes_list,
		barcodes_fastq = barcodes_fastq,
		reads_fastq = reads_fastq
	}

	call bwaAlignmentTask {
		input:
		file_glob=splitBarcodesTask.reads,
		ref_idxs=ref_idxs,
		ref_fasta=ref_fasta
	}

	call indelAccuracyTask {
		input:
		bams=bwaAlignmentTask.bams,
		idxs=bwaAlignmentTask.idxs,
		indel_cut_interval=indel_cut_interval,
		indel_cut_site=indel_cut_site,
		indel_genes=indel_genes,
		ref_idxs=ref_idxs,
		ref_fasta=ref_fasta
	}

	call novelIndelTask {
		input:
		bams=bwaAlignmentTask.bams,
		idxs=bwaAlignmentTask.idxs,
		indel_cut_interval=indel_cut_interval,
		indel_genes=indel_genes,
		neg_cntrl=neg_cntrl
	}
}

task splitBarcodesTask {
	File barcodes_list
	File barcodes_fastq
	File reads_fastq

	command <<<
	Rscript /usr/src/app/src/SplitBarcodes.R . ${barcodes_fastq} ${reads_fastq} ${barcodes_list}
	>>>

	output {
		Array[File] reads = glob("Reads/*.fastq")
		File mappings = "barcode_mapping.csv"
		File unmapped_reads = "unknown_reads.fastq"
		File unmapped_barcodes = "unknown_barcodes.fastq"
	}

	runtime {
		docker: "mburger/crispr-seq"
	}

}

task bwaAlignmentTask {
	Array[File] file_glob
	File ref_idxs
	Array[Array[File]] refFiles = read_tsv(ref_idxs)
	File ref_fasta

	command <<<
		mkdir Reads
		cd Reads
		/usr/src/app/src/bwaAlignmentTask.sh ${sep=',' file_glob} ${ref_fasta}
		cd ..  
	>>>

	output {
		Array[File] bams = glob("Reads/*.bam")
		Array[File] idxs = glob("Reads/*.bai")
	}

	runtime {
		docker: "mburger/crispr-seq"
	}
}

task novelIndelTask {
	Array[File] bams
	Array[File] idxs
	File indel_cut_interval
	File indel_genes
	File neg_cntrl

	command <<<
		mkdir Reads
		/usr/src/app/src/mvFiles.sh ${sep=',' bams} ./Reads
		/usr/src/app/src/mvFiles.sh ${sep=',' idxs} ./Reads
		/usr/src/app/src/novelIndelTask.sh /usr/src/app/src . ${indel_cut_interval} ${indel_genes} ${neg_cntrl}
		tar -cf VariantCalls.tar VariantCalls/
		tar -cf IndelQuant.tar indel_quant/
	>>>

	output {
		File VC = "VariantCalls.tar"
		File IQ = "IndelQuant.tar"
	}

	runtime {
		docker: "mburger/crispr-seq"
	}

}

task indelAccuracyTask {
	Array[File] bams
	Array[File] idxs
	File indel_cut_site
	File indel_cut_interval
	File indel_genes

	File ref_idxs
	Array[Array[File]] refFiles = read_tsv(ref_idxs)
	File ref_fasta

	command <<<
		mkdir Reads
		/usr/src/app/src/mvFiles.sh ${sep=',' bams} ./Reads
		/usr/src/app/src/mvFiles.sh ${sep=',' idxs} ./Reads
		/usr/src/app/src/indelAccuracyTask.sh . /usr/src/app/src ${indel_cut_site} ${indel_genes} ${ref_fasta} ${indel_cut_interval}
		tar -cf Power.tar Power/
	>>>

	output {
		File QA = "Power.tar"
	}

	runtime {
		docker: "mburger/crispr-seq"
	}

}


workflow crisprSeqWorkflow {
	File barcodes_list
	File barcodes_fastq
	File reads_fastq

	File ref_idxs
	File ref_fasta

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
		#ref_amb=ref_amb,
		ref_fasta=ref_fasta
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


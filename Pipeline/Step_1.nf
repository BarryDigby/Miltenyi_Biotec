#!/usr/bin/env nextflow

files = Channel.fromPath("/data/bdigby/Projects/Galway_Genomics/Files/*.fastq.gz")

process fastqc {
        publishDir "fastqc/Pre-Trim/", mode: 'copy'

        input:
        file reads from files

        output:
        file "*_fastqc.{zip,html}" into fastqc_results

        script:
        """
        fastqc -q ${reads}
        """
}

process multiqc {
        publishDir "fastqc/Pre-Trim/", mode: 'copy'

        input:
        file ('*') from fastqc_results.collect().ifEmpty([])

        output:
        file "multiqc_report.html" into multiqc_report
        file "multiqc_data"

        script:
        """
        multiqc .
        """
}

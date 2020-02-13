#!/usr/bin/env nextflow

params.reads = "Files/*R{1,2}_001.fastq.gz"
Channel
        .fromFilePairs( params.reads )
        .set { input_fq }

adapter_file = "/repo/bbmap/resources/adapters.fa"

process bbduk {
        publishDir "trimmed_reads/", mode: 'copy'

        input:
        set val(pair_id), file(reads) from input_fq
        file adapter_file

        output:
        set val(pair_id), file("*trimmed.fastq.gz") into fastqc_input
        file "${pair_id}.stats.txt"

        script:
        """
        echo "read pair ${pair_id}"
        echo "reads[0] - ${reads[0]}"
        echo "reads[0].baseName - ${reads[0].baseName}"
        bbduk \
        -Xmx4g \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1=${pair_id}_r1.trimmed.fastq.gz \
        out2=${pair_id}_r2.trimmed.fastq.gz \
        literal=AGATCGGAAGAG \
        minlen=30 \
        ktrim=r \
        k=12 \
        qtrim=rl \
        trimq=20 \
        stats=${pair_id}.stats.txt
        """
}

process fastqc {
        publishDir "fastqc/Post-Trim/", mode: 'copy'

        input:
        set val(pair_id), file(trimmed_reads) from fastqc_input

        output:
        file "*.{zip,html}" into fastqc_output

        script:
        """
        fastqc -t 8 -q ${trimmed_reads}
        """
}

process multiqc {
        publishDir "fastqc/Post-Trim/", mode:'copy'

        input:
        file("*") from fastqc_output.collect().ifEmpty([])

        output:
        file "multiqc_report.html" into multiqc_report
        file "multiqc_data"

        script:
        """
        multiqc .
        """
}

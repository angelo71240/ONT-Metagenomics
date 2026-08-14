process FILTLONG {

    tag "$sample"

    publishDir "${params.outdir}/filtlong", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_filtered.fastq.gz"), emit: filtered_reads

    script:
    """
    filtlong \
        --min_length 500 \
        --keep_percent 90 \
        ${reads} | gzip > ${sample}_filtered.fastq.gz
    """
}

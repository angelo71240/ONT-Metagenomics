process KRAKEN2 {

    tag "$sample"

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample),
          path("${sample}.kraken.out"),
          path("${sample}.kraken.report")

    script:
    """
    kraken2 \
        --db ${params.kraken_db} \
        --gzip-compressed \
        --output ${sample}.kraken.out \
        --report ${sample}.kraken.report \
        ${reads}
    """
}

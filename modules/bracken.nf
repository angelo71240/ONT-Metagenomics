process BRACKEN {

    tag "$sample"

    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(sample), path(kraken_report)

    output:
    tuple val(sample),
          path("${sample}.bracken.tsv")

    script:
    """
    bracken \
        -d ${params.kraken_db} \
        -i ${kraken_report} \
        -o ${sample}.bracken.tsv \
        -r 300 \
        -l S
    """
}

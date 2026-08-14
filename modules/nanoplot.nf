process NANOPLOT {

    tag "$sample"

    publishDir "${params.outdir}/nanoplot/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_nanoplot"), emit: reports

    script:
    """
    mkdir -p ${sample}_nanoplot

    NanoPlot \
        --fastq ${reads} \
        --outdir ${sample}_nanoplot \
        --prefix ${sample} \
        --threads ${task.cpus}
    """
}

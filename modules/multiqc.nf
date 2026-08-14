process MULTIQC {

    tag "MultiQC"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . \
        --outdir . \
        --filename multiqc_report.html
    """
}

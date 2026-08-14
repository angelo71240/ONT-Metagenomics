process DIFFERENTIAL_ABUNDANCE {

    tag "Abundance comparison"

    publishDir "${params.outdir}/differential", mode: 'copy'

    input:
    path abundance_matrix

    output:
    path "differential_results.tsv"
    path "abundance_barplot.png"
    path "abundance_heatmap.png"

    script:
    """
    python3 ${projectDir}/scripts/differential_abundance.py
    """
}

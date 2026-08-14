process ABUNDANCE_MATRIX {

    tag "Abundance matrix"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path bracken_files

    output:
    path "abundance_matrix.tsv"

    script:
    """
    python3 ${projectDir}/scripts/make_abundance_matrix.py \
        ${bracken_files.join(' ')}
    """
}

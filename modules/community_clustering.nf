process COMMUNITY_CLUSTERING {

    tag "Community clustering"

    publishDir "${params.outdir}/community_clustering", mode: 'copy'

    input:
    path abundance_matrix

    output:
    path "bray_curtis_distance.tsv"
    path "pcoa_coordinates.tsv"
    path "pcoa_bray_curtis.png"
    path "hierarchical_clustering.png"
    path "pca_coordinates.tsv"
    path "pca.png"
    path "community_heatmap.png"

    script:
    """
    python3 ${projectDir}/scripts/community_clustering.py
    """
}

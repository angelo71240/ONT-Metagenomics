nextflow.enable.dsl=2

params.reads     = "data/*.fastq.gz"
params.outdir    = "results"
params.kraken_db = "/mnt/scratch/DR/metagenome-pipeline/assets/databases/kraken2"


include { NANOPLOT as NANOPLOT_PRE }  from './modules/nanoplot.nf'
include { FILTLONG }                  from './modules/filtlong.nf'
include { NANOPLOT as NANOPLOT_POST } from './modules/nanoplot.nf'
include { MULTIQC }                   from './modules/multiqc.nf'

include { KRAKEN2 }                   from './modules/kraken2.nf'
include { BRACKEN }                   from './modules/bracken.nf'

include { ABUNDANCE_MATRIX }          from './modules/abundance_matrix.nf'
include { DIFFERENTIAL_ABUNDANCE }    from './modules/differential_abundance.nf'
include { COMMUNITY_CLUSTERING }      from './modules/community_clustering.nf'


workflow {

    /*
     * STEP 1 - INPUT
     * Load raw ONT FASTQ files
     */

    reads_ch = Channel
        .fromPath(params.reads)
        .map { reads ->
            tuple(reads.simpleName, reads)
        }


    /*
     * STEP 2 - NANOPLOT PRE-FILTER QC
     * QC of raw ONT reads
     */

    NANOPLOT_PRE(reads_ch)


    /*
     * STEP 3 - FILTLONG
     * Filter short and/or low-quality ONT reads
     */

    FILTLONG(reads_ch)


    /*
     * STEP 4 - NANOPLOT POST-FILTER QC
     * QC of filtered ONT reads
     */

    nanoplot_post_input = FILTLONG.out.filtered_reads

    NANOPLOT_POST(nanoplot_post_input)


    /*
     * STEP 5 - MULTIQC
     * Summarize NanoPlot QC results
     */

    multiqc_input = NANOPLOT_PRE.out.reports
        .mix(NANOPLOT_POST.out.reports)
        .collect()

    MULTIQC(multiqc_input)


    /*
     * STEP 6 - KRAKEN2
     * Taxonomic classification of filtered ONT reads
     */

    KRAKEN2(FILTLONG.out.filtered_reads)


    /*
     * STEP 7 - BRACKEN
     * Taxonomic abundance estimation
     */

    bracken_input = KRAKEN2.out.map { sample, kraken_out, kraken_report ->
        tuple(sample, kraken_report)
    }

    BRACKEN(bracken_input)


    /*
     * STEP 8 - ABUNDANCE MATRIX
     * Combine Bracken results into one matrix
     */

    bracken_files = BRACKEN.out
        .map { sample, bracken_file -> bracken_file }
        .collect()

    ABUNDANCE_MATRIX(bracken_files)


    /*
     * STEP 9 - DIFFERENTIAL ABUNDANCE (project Step 5)
     * Calculate abundance differences and generate plots
     */

    DIFFERENTIAL_ABUNDANCE(ABUNDANCE_MATRIX.out)


    /*
     * STEP 10 - COMMUNITY CLUSTERING (project Step 6)
     * Bray-Curtis, PCoA, hierarchical clustering, PCA, heatmap
     */

    COMMUNITY_CLUSTERING(ABUNDANCE_MATRIX.out)
}

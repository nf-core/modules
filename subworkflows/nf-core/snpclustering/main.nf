include { BEAGLE5_BEAGLE                } from '../../../modules/nf-core/beagle5/beagle'
include { PLINK2_VCF                    } from '../../../modules/nf-core/plink2/vcf'
include { PLINK2_PCA                    } from '../../../modules/nf-core/plink2/pca'
include { GAWK as GAWK_EIGENVEC_TO_TSV  } from '../../../modules/nf-core/gawk'
include { CUSTOM_PCACLUSTERING          } from '../../../modules/nf-core/custom/pcaclustering'
include { CUSTOM_CLUSTERMETRICS         } from '../../../modules/nf-core/custom/clustermetrics'
include { CUSTOM_CLUSTERVISUALIZATION   } from '../../../modules/nf-core/custom/clustervisualization'

workflow SNPCLUSTERING {
    take:
    vcf_ch
    refpanel_ch
    genmap_ch
    region
    n_pcs
    use_approx
    algorithm
    n_clusters
    dbscan_eps
    dbscan_min_samples

    main:
    // Bundled awk program (see ./awk/eigenvec_to_tsv.awk); shipped alongside
    // the subworkflow and passed to GAWK_EIGENVEC_TO_TSV as `program_file` so
    // the awk source stays readable as awk rather than as an escaped Groovy
    // string.
    ch_eigenvec_to_tsv_awk = file("${moduleDir}/awk/eigenvec_to_tsv.awk", checkIfExists: true)

    /*
     * Build BEAGLE input tuple:
     * tuple val(meta), path(vcf), path(vcf_index), path(refpanel), path(refpanel_index),
     *       path(genmap), path(exclsamples), path(exclmarkers), val(region)
     */
    ch_beagle_input = vcf_ch.map { meta, vcf, vcf_index ->
        tuple(
            meta,
            vcf,
            vcf_index,
            [],
            [],
            [],
            [],
            [],
            region
        )
    }

    BEAGLE5_BEAGLE(ch_beagle_input)

    /*
     * Convert imputed VCF to PLINK2 pfiles
     */
    PLINK2_VCF(BEAGLE5_BEAGLE.out.vcf)

    /*
     * PLINK2_PCA expects:
     * tuple val(meta), val(npcs), val(use_approx), path(pgen), path(psam), path(pvar)
     */
    ch_plink_pca_input = PLINK2_VCF.out.pgen
        .join(PLINK2_VCF.out.pvar)
        .join(PLINK2_VCF.out.psam)
        .map { meta, pgen, pvar, psam ->
            tuple(meta, npcs, use_approx, pgen, psam, pvar)
        }

    PLINK2_PCA(ch_plink_pca_input)

    /*
     * Convert .eigenvec to a tab-separated file for downstream clustering
     * modules: strips the leading '#' from the header and renames the IID
     * column to sample_id. Output keeps its default naming (<meta.id>.eigenvec)
     * — only the tab-separated content matters downstream, not the extension.
     */
    GAWK_EIGENVEC_TO_TSV(
        PLINK2_PCA.out.evecfile,
        ch_eigenvec_to_tsv_awk,
        false
    )

    /*
     * PCA clustering
     * CUSTOM_PCACLUSTERING(tsv, algorithm, n_clusters, dbscan_eps, dbscan_min_samples)
     */
    CUSTOM_PCACLUSTERING(
        GAWK_EIGENVEC_TO_TSV.out.output,
        algorithm,
        n_clusters,
        dbscan_eps,
        dbscan_min_samples
    )

    /*
     * Metrics and visualization both expect one tuple input channel
     * built from: meta + tsv + cluster assignments
     */
    ch_cluster_analysis_input = GAWK_EIGENVEC_TO_TSV.out.output
        .join(CUSTOM_PCACLUSTERING.out.clusters)
        .map { meta, tsv, clusters ->
            tuple(meta, tsv, clusters)
        }

    CUSTOM_CLUSTERMETRICS(ch_cluster_analysis_input)
    CUSTOM_CLUSTERVISUALIZATION(ch_cluster_analysis_input)

    emit:
    imputed_vcf   = BEAGLE5_BEAGLE.out.vcf
    beagle_log    = BEAGLE5_BEAGLE.out.log
    pgen          = PLINK2_VCF.out.pgen
    pvar          = PLINK2_VCF.out.pvar
    psam          = PLINK2_VCF.out.psam
    evecfile      = PLINK2_PCA.out.evecfile
    evfile        = PLINK2_PCA.out.evfile
    pca_log       = PLINK2_PCA.out.logfile
    tsv           = GAWK_EIGENVEC_TO_TSV.out.output
    clusters      = CUSTOM_PCACLUSTERING.out.clusters
    cluster_info  = CUSTOM_PCACLUSTERING.out.info
    metrics       = CUSTOM_CLUSTERMETRICS.out.metrics
    k_sweep       = CUSTOM_CLUSTERMETRICS.out.k_sweep
    selected      = CUSTOM_CLUSTERMETRICS.out.selected
    metric_plots  = CUSTOM_CLUSTERMETRICS.out.plots
    umap_tsv      = CUSTOM_CLUSTERVISUALIZATION.out.umap_tsv
    tsne_tsv      = CUSTOM_CLUSTERVISUALIZATION.out.tsne_tsv
    umap_png      = CUSTOM_CLUSTERVISUALIZATION.out.umap_png
    tsne_png      = CUSTOM_CLUSTERVISUALIZATION.out.tsne_png
}

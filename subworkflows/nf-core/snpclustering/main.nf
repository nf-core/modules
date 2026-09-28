include { BEAGLE5_BEAGLE              } from '../../../modules/nf-core/beagle5/beagle'
include { PLINK2_VCF                  } from '../../../modules/nf-core/plink2/vcf'
include { PLINK2_PCA                  } from '../../../modules/nf-core/plink2/pca'
include { CUSTOM_PCACLUSTERING        } from '../../../modules/nf-core/custom/pcaclustering'
include { CUSTOM_CLUSTERMETRICS       } from '../../../modules/nf-core/custom/clustermetrics'
include { CUSTOM_CLUSTERVISUALIZATION } from '../../../modules/nf-core/custom/clustervisualization'

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
    PLINK2_VCF(BEAGLE5_BEAGLE.out.vcf)

    ch_plink_pca_input = PLINK2_VCF.out.pgen
        .join(PLINK2_VCF.out.pvar)
        .join(PLINK2_VCF.out.psam)
        .map { meta, pgen, pvar, psam ->
            tuple(meta, n_pcs, use_approx, pgen, psam, pvar)
        }

    PLINK2_PCA(ch_plink_pca_input)

    CUSTOM_PCACLUSTERING(
        PLINK2_PCA.out.evecfile,
        algorithm,
        n_clusters,
        dbscan_eps,
        dbscan_min_samples
    )

    ch_cluster_analysis_input = PLINK2_PCA.out.evecfile
        .join(CUSTOM_PCACLUSTERING.out.clusters)
        .map { meta, eigenvec, clusters ->
            tuple(meta, eigenvec, clusters)
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

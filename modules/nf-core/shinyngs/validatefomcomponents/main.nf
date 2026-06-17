process SHINYNGS_VALIDATEFOMCOMPONENTS {
    tag "$sample"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d782b4f11adf8f3cad6af74ea585468decd873a171da1dae0e4a24a82bb29020/data' :
        'community.wave.seqera.io/library/r-shinyngs:2.4.0--709fc6932be670a5' }"

    input:
    tuple val(meta),  path(sample), path(assay_files)
    tuple val(meta2), path(feature_meta)
    tuple val(meta3), path(contrasts)

    output:
    tuple val(meta), path("*/*.sample_metadata.tsv")   , emit: sample_meta
    tuple val(meta), path("*/*.feature_metadata.tsv")  , emit: feature_meta, optional: true
    tuple val(meta), path("*/*.assay.tsv")             , emit: assays
    tuple val(meta), path("*/*.contrasts_file.tsv")    , emit: contrasts
    tuple val("${task.process}"), val('shinyngs'), eval('Rscript -e "library(shinyngs); cat(as.character(packageVersion(\'shinyngs\')))"'), emit: versions_shinyngs, topic: versions


    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/validate_fom_components.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def feature = feature_meta ? "--feature_metadata '$feature_meta'" : ''

    """
    validate_fom_components.R \\
        --sample_metadata "$sample" \\
        $feature \\
        --assay_files "${assay_files.join(',')}" \\
        --contrasts_file "$contrasts" \\
        --output_directory "$prefix" \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p $prefix

    printf 'sample\\tcondition\\nsample1\\tcondition1\\nsample2\\tcondition2\\n' \\
        > $prefix/${prefix}.sample_metadata.tsv

    printf 'gene_id\\tgene_name\\ngene1\\tname1\\n' \\
        > $prefix/${prefix}.feature_metadata.tsv

    printf 'gene_id\\tsample1\\tsample2\\ngene1\\t1\\t2\\n' \\
        > $prefix/${prefix}.assay.tsv

    printf 'id\\tvariable\\treference\\ttarget\\tblocking\\tformula\\tmake_contrasts_str\\ncontrast1\\tcondition\\tcondition1\\tcondition2\\t\\t\\t\\n' \\
        > $prefix/${prefix}.contrasts_file.tsv
    """
}

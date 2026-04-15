process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2df03eec9c0805810e0ef6caec4347d7c6545eece61e941018945502fafc9b6/data'
        : 'community.wave.seqera.io/library/gtdbtk_python:cee0379cf1ca2968'}"

    input:
    tuple val(meta)   , path("bins/*")
    tuple val(db_name), path(db)
    val use_pplacer_scratch_dir

    output:
    tuple val(meta), path("${prefix}")                               , emit: gtdb_outdir
    tuple val(meta), path("${prefix}/classify/*.summary.tsv")        , emit: summary
    tuple val(meta), path("${prefix}/classify/*.classify.tree")      , emit: tree       , optional: true
    tuple val(meta), path("${prefix}/identify/*.markers_summary.tsv"), emit: markers    , optional: true
    tuple val(meta), path("${prefix}/align/*.msa.fasta.gz")          , emit: msa        , optional: true
    tuple val(meta), path("${prefix}/align/*.user_msa.fasta.gz")     , emit: user_msa   , optional: true
    tuple val(meta), path("${prefix}/align/*.filtered.tsv")          , emit: filtered   , optional: true
    tuple val(meta), path("${prefix}/identify/*.failed_genomes.tsv") , emit: failed     , optional: true
    tuple val(meta), path("${prefix}/${prefix}.log")                 , emit: log
    tuple val(meta), path("${prefix}/${prefix}.warnings.log")        , emit: warnings
    tuple val("${task.process}"), val('gtdbtk'), eval("gtdbtk --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1") , topic: versions, emit: versions_gtdbtk
    tuple val("${task.process}"), val('gtdb_db'), eval('grep VERSION_DATA $GTDBTK_DATA_PATH/metadata/metadata.txt | sed "s/VERSION_DATA=//"'), topic: versions, emit: versions_gtdbtk_db

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    prefix              = task.ext.prefix ?: "${meta.id}"
    def pplacer_scratch = use_pplacer_scratch_dir ? "--scratch_dir pplacer_tmp" : ""
    """
    export GTDBTK_DATA_PATH="\$(find -L ${db} -name 'metadata' -type d -exec dirname {} \\;)"

    if [ "${pplacer_scratch}" != "" ] ; then
        mkdir pplacer_tmp
    fi

    gtdbtk classify_wf \\
        ${args} \\
        --genome_dir bins \\
        --prefix "${prefix}" \\
        --out_dir ${prefix} \\
        --cpus ${task.cpus} \\
        ${pplacer_scratch}

    mv ${prefix}/gtdbtk.log "${prefix}/${prefix}.log"
    mv ${prefix}/gtdbtk.warnings.log "${prefix}/${prefix}.warnings.log"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    mkdir ${prefix}/identify
    mkdir ${prefix}/classify
    mkdir ${prefix}/align

    touch ${prefix}/classify/${prefix}.ar53.summary.tsv
    touch ${prefix}/classify/${prefix}.bac120.summary.tsv
    touch ${prefix}/classify/${prefix}.ar53.classify.tree
    touch ${prefix}/classify/${prefix}.bac120.classify.tree

    touch ${prefix}/identify/${prefix}.ar53.markers_summary.tsv
    touch ${prefix}/identify/${prefix}.bac120.markers_summary.tsv

    echo "" | gzip > ${prefix}/align/${prefix}.ar53.msa.fasta.gz
    echo "" | gzip > ${prefix}/align/${prefix}.bac120.user_msa.fasta.gz
    touch ${prefix}/align/${prefix}.ar53.filtered.tsv
    touch ${prefix}/align/${prefix}.bac120.filtered.tsv

    touch ${prefix}/${prefix}.log
    touch ${prefix}/${prefix}.warnings.log
    touch ${prefix}/${prefix}.failed_genomes.tsv
    """
}

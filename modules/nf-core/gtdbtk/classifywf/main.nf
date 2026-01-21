process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.5.2--pyh1f0d9b5_0':
        'biocontainers/gtdbtk:2.5.2--pyh1f0d9b5_0' }"

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
    path ("versions.yml")                                            , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version 2>/dev/null) | sed "s/gtdbtk: version //; s/ Copyright.*//")
        gtdb_db: \$(grep VERSION_DATA \$GTDBTK_DATA_PATH/metadata/metadata.txt | sed "s/VERSION_DATA=//")
    END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version 2>/dev/null) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}

process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.4.1--pyhdfd78af_1':
        'biocontainers/gtdbtk:2.4.1--pyhdfd78af_1' }"

    input:
    tuple val(meta)   , path("bins/*")
    tuple val(db_name), path(db)
    val use_pplacer_scratch_dir
    path mash_db

    output:
    tuple val(meta), path("*")                             , emit: gtdb_files
    tuple val(meta), path("classify/*.summary.tsv")        , emit: summary
    tuple val(meta), path("classify/*.classify.tree")      , emit: tree       , optional: true
    tuple val(meta), path("identify/*.markers_summary.tsv"), emit: markers    , optional: true
    tuple val(meta), path("align/*.msa.fasta.gz")          , emit: msa        , optional: true
    tuple val(meta), path("align/*.user_msa.fasta.gz")     , emit: user_msa   , optional: true
    tuple val(meta), path("align/*.filtered.tsv")          , emit: filtered   , optional: true
    tuple val(meta), path("identify/*.failed_genomes.tsv") , emit: failed     , optional: true
    tuple val(meta), path("gtdbtk.${prefix}.log")          , emit: log
    tuple val(meta), path("gtdbtk.${prefix}.warnings.log") , emit: warnings
    path ("versions.yml")                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    prefix              = task.ext.prefix ?: "${meta.id}"
    def pplacer_scratch = use_pplacer_scratch_dir ? "--scratch_dir pplacer_tmp" : ""
    def mash_mode       = mash_db ? "--mash_db ${mash_db}" : "--skip_ani_screen"
    """
    export GTDBTK_DATA_PATH="\$(find -L ${db} -name 'metadata' -type d -exec dirname {} \\;)"

    if [ "${pplacer_scratch}" != "" ] ; then
        mkdir pplacer_tmp
    fi

    gtdbtk classify_wf \\
        ${args} \\
        --genome_dir bins \\
        --prefix "gtdbtk.${prefix}" \\
        --out_dir "\${PWD}" \\
        --cpus ${task.cpus} \\
        ${mash_mode} \\
        ${pplacer_scratch}

    mv gtdbtk.log "gtdbtk.${prefix}.log"
    mv gtdbtk.warnings.log "gtdbtk.${prefix}.warnings.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version 2>/dev/null) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir identify
    mkdir classify
    mkdir align

    touch classify/gtdbtk.${prefix}.ar53.summary.tsv
    touch classify/gtdbtk.${prefix}.bac120.summary.tsv
    ln -s classify/gtdbtk.${prefix}.ar53.summary.tsv gtdbtk.${prefix}.ar53.summary.tsv
    ln -s classify/gtdbtk.${prefix}.bac120.summary.tsv gtdbtk.${prefix}.bac120.summary.tsv
    touch classify/gtdbtk.${prefix}.ar53.classify.tree
    touch classify/gtdbtk.${prefix}.bac120.classify.tree

    touch identify/gtdbtk.${prefix}.ar53.markers_summary.tsv
    touch identify/gtdbtk.${prefix}.bac120.markers_summary.tsv

    echo "" | gzip > align/gtdbtk.${prefix}.ar53.msa.fasta.gz
    echo "" | gzip > align/gtdbtk.${prefix}.bac120.user_msa.fasta.gz
    touch align/gtdbtk.${prefix}.ar53.filtered.tsv
    touch align/gtdbtk.${prefix}.bac120.filtered.tsv

    touch gtdbtk.${prefix}.log
    touch gtdbtk.${prefix}.warnings.log
    touch gtdbtk.${prefix}.failed_genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version 2>/dev/null) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}

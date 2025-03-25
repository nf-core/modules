process INFERNAL_CMSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h7b50bb2_4' :
        'biocontainers/infernal:1.1.5--pl5321h7b50bb2_4' }"

    input:
    tuple val(meta), path(cmfile), path(seqdb)
    val(write_align)
    val(write_target)

    output:
    tuple val(meta), path('*.txt.gz'), emit: output
    tuple val(meta), path('*.sto.gz'), emit: alignments    , optional: true
    tuple val(meta), path('*.tbl.gz'), emit: target_summary, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    output          = "${prefix}.txt"
    alignment       = write_align     ? "-A ${prefix}.sto" : ''
    target_summary  = write_target    ? "--tblout ${prefix}.tbl" : ''
    def seqdb_input = seqdb.toString() - ~/\.gz$/
    def gunzip      = seqdb.getExtension() == "gz" ? "gunzip -c ${seqdb} > ${seqdb_input}" : ""
    def cleanup     = seqdb.getExtension() == "gz" ? "rm ${seqdb_input}" : ""
    """
    ${gunzip}

    cmsearch \\
        $args \\
        --cpu ${task.cpus} \\
        -o ${output} \\
        ${alignment} \\
        ${target_summary} \\
        ${cmfile} \\
        ${seqdb_input}

    gzip --no-name *.txt \\
        ${write_align ? '*.sto' : ''} \\
        ${write_target ? '*.tbl' : ''}

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmsearch -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.txt"
    ${write_align ? "touch ${prefix}.sto" : ''} \\
    ${write_target ? "touch ${prefix}.tbl" : ''}

    gzip --no-name *.txt \\
        ${write_align ? '*.sto' : ''} \\
        ${write_target ? '*.tbl' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmsearch -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//')
    END_VERSIONS
    """
}

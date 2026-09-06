process INFERNAL_CMSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56cfd2ffd3ea93fb3611175e7383ca253108ad32e868bc08067bf9197da408af/data' :
        'community.wave.seqera.io/library/infernal:1.1.5--8877beae3740ff72' }"

    input:
    tuple val(meta), path(cmfile), path(seqdb)
    val(write_align)
    val(write_target)

    output:
    tuple val(meta), path('*.txt.gz'), emit: output
    tuple val(meta), path('*.sto.gz'), emit: alignments    , optional: true
    tuple val(meta), path('*.tbl.gz'), emit: target_summary, optional: true
    tuple val("${task.process}"), val('infernal'), eval("cmsearch -h | sed '2!d;s/.*INFERNAL //;s/ .*//'"), emit: versions_infernal, topic: versions
    tuple val("${task.process}"), val('gzip'), eval("gzip --version |& sed '1!d;s/gzip //'"), topic: versions, emit: versions_gzip

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.txt.gz
    ${write_align  ? "echo '' | gzip > ${prefix}.sto.gz" : ''}
    ${write_target ? "echo '' | gzip > ${prefix}.tbl.gz" : ''}
    """
}

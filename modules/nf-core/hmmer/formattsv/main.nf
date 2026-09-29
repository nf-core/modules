process HMMER_FORMATTSV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/0600fdb9287986d3d684b95d6e94f597c41b24a746608532d7c0836bd7b44c15/data'
:         'community.wave.seqera.io/library/gawk_gzip:34c193f8cae8819b' }"

    input:
    tuple val(meta), val(labels), path(files, arity: '1..*')
    val format // 'tblout' or 'domtblout' -- selects which HMMER column layout to parse

    output:
    tuple val(meta), path("*.tsv.gz"),                                          emit: tsv
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    assert labels.size() == files.size() : "HMMER_FORMATTSV: labels and files must be the same length (got ${labels.size()} labels, ${files.size()} files)"
    assert format in ['tblout', 'domtblout'] : "HMMER_FORMATTSV: format must be 'tblout' or 'domtblout', got '${format}'"
    if (format == 'tblout') {
        template 'format_tblout.sh'
    } else {
        template 'format_domtblout.sh'
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    assert format in ['tblout', 'domtblout'] : "HMMER_FORMATTSV: format must be 'tblout' or 'domtblout', got '${format}'"
    """
    echo "" | gzip > ${prefix}.${format}.tsv.gz
    """
}

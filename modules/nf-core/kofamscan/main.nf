// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process KOFAMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/98/98c22912f6c82843ff35f6fd072f9329e821da1a897cdf44cc2f812d3ca999e8/data'
        : 'community.wave.seqera.io/library/kofamscan:1.3.0--0a5b7fcd8f6e6678'}"

    input:
    tuple val(meta), path(fasta)
    path profiles
    path ko_list

    output:
    tuple val(meta), path("${prefix}.txt.gz"), optional: true, emit: txt
    tuple val(meta), path("${prefix}.tsv.gz"), optional: true, emit: tsv
    tuple val("${task.process}"), val('kofamscan'), eval("exec_annotation --version 2>&1 | sed 's/exec_annotation //;'"), topic: versions, emit: versions_kofamscan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--format detail-tsv") ? "tsv" :
                    "txt"
    def is_compressed = fasta.name.endsWith('.gz')
    def fasta_name = is_compressed ? fasta.name - '.gz' : fasta.name
    def gunzip = is_compressed ? "gunzip -c ${fasta} > ${fasta_name}" : ''
    def cleanup = is_compressed ? "rm -f ${fasta_name}" : ''
    """
    $gunzip

    exec_annotation \\
        -p $profiles \\
        -k $ko_list \\
        $args \\
        --cpu $task.cpus \\
        -o ${prefix}.${extension} \\
        ${fasta_name}

    gzip -n ${prefix}.${extension}

    $cleanup
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--format detail-tsv") ? "tsv" :
                    "txt"
    """
    touch ${prefix}.${extension}
    gzip -n ${prefix}.${extension}
    """
}

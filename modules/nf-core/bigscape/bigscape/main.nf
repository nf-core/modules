process BIGSCAPE_BIGSCAPE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ac/acb20daff43670baa2f5599e97ef0d96f3eeeaaaeb31e7dfa7a23066c908b44a/data' :
        'community.wave.seqera.io/library/bigscape:2.0.3--aed310da2cd1ca1d' }"

    input:
    tuple val(meta), path(gbk_dir)
    path pfam_hmm

    output:
    tuple val(meta), path("${prefix}"), emit: output
    tuple val("${task.process}"), val('bigscape'), eval("bigscape --version | sed 's/BiG-SCAPE //'"), topic: versions, emit: versions_bigscape

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def pfam_file = pfam_hmm.name.endsWith('.gz') ? pfam_hmm.baseName : pfam_hmm
    """
    if [[ "${pfam_hmm}" == *.gz ]]; then
        gzip -dc ${pfam_hmm} > ${pfam_file}
    fi

    bigscape cluster \\
        --gbk-dir ${gbk_dir} \\
        --pfam-path ${pfam_file} \\
        --output-dir ${prefix} \\
        --cores ${task.cpus} \\
        ${args}

    if [[ "${pfam_hmm}" == *.gz ]]; then
        rm -f ${pfam_file}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/network_files
    touch ${prefix}/network_files/mix.network
    mkdir -p ${prefix}/html_content
    touch ${prefix}/html_content/index.html
    """
}

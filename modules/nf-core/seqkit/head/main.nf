process SEQKIT_HEAD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(fastqs), val(seq_count)

    output:
    tuple val(meta), path("${prefix}_subset_*"), emit: subset
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    for f in ${fastqs.join(' ')}
    do
        seqkit head \\
            ${args} \\
            --threads ${task.cpus} \\
            -n ${seq_count} \\
            -o "${prefix}_subset_\$(basename \$f)" \\
            \$f
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in ${fastqs.join(' ')}
    do
       echo "" | gzip > "${prefix}_subset_\$(basename \$f)"
    done
    """
}

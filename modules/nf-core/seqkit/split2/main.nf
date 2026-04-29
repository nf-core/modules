process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data' :
        'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}/*"), emit: reads
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        seqkit \\
            split2 \\
            $args \\
            --threads $task.cpus \\
            $reads \\
            --out-dir ${prefix}
        """
    } else {
        """
        seqkit \\
            split2 \\
            $args \\
            --threads $task.cpus \\
            --read1 ${reads[0]} \\
            --read2 ${reads[1]} \\
            --out-dir ${prefix}
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        mkdir -p ${prefix}
        echo "" | gzip > ${prefix}/${reads[0]}
        """
    } else {
        """
        mkdir -p ${prefix}
        echo "" | gzip > ${prefix}/${reads[0]}
        echo "" | gzip > ${prefix}/${reads[1]}
        """
    }
}

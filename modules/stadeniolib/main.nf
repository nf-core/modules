process STADENIOLIB {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::staden_io_lib=1.14.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staden_io_lib:1.14.14--h9dace67_0' :
        'quay.io/biocontainers/staden_io_lib:1.14.14--h9dace67_0' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.cram") ,emit: cram
    path "versions.yml"             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    scramble \
        -I bam \
        -O cram \
        -r ${fasta} \
        -t $task.cpus \
        ${bam} \
        ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stadeniolib: \$(echo \$(scramble -h 2>&1) | sed 's/^.*version //; s/Author.*\$//' )
    END_VERSIONS
    """
}

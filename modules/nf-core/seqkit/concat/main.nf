process SEQKIT_CONCAT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(input, stageAs: 'in/*')

    output:
    tuple val(meta), path("*.{fasta,fastq,fa,fq,fas,fna,faa}"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input instanceof List ? input[0].getExtension() : input.getExtension()
    """
    seqkit \\
        concat \\
        --threads ${task.cpus} \\
        ${args} \\
        in/* > ${prefix}.${file_type}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    """
}

process HOLODECK_METHYLATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fec611479d4b632095f4173ab61e9363051bd5645386a7873ce378846461cae9/data'
        : 'community.wave.seqera.io/library/holodeck:0.3.0--f29a8a9e2f667299'}"

    input:
    tuple val(meta), path(fasta), path(fai), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('holodeck'), eval("holodeck --version | sed 's/^holodeck //'"), emit: versions_holodeck, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_arg = vcf ? "--vcf ${vcf}" : ""
    """
    holodeck \\
        methylate \\
        --reference ${fasta} \\
        ${vcf_arg} \\
        --output ${prefix}.vcf.gz \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}

process TRGT_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b94772fdcba7aa8027797b0d092a77b272cebc8d566fec3c960a1470e20223cf/data':
        'community.wave.seqera.io/library/trgt:5.1.0--66e3b26326e566e7' }"

    input:
    tuple val(meta) , path(bam), path(bai), val(karyotype)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(repeats)

    output:
    tuple val(meta), path("*.vcf.gz")      , emit: vcf
    tuple val(meta), path("*.spanning.bam"), emit: bam, optional: true
    tuple val("${task.process}"), val('trgt'), eval("trgt --version | sed 's/.* //g'"), emit: versions_trgt, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def karyo = karyotype ? "--karyotype ${karyotype}" : ""
    """
    trgt genotype \\
        $args \\
        --genome ${fasta} \\
        --reads ${bam} \\
        --repeats ${repeats} \\
        ${karyo} \\
        --threads ${task.cpus} \\
        --output-prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.spanning.bam
    echo "" | gzip > ${prefix}.vcf.gz
    """
}

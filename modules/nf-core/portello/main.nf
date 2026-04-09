process PORTELLO {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b7a7ecbc4e475e2b0a339facbda324068cab503a58261d6a235169ce04bf25d8/data'
        : 'community.wave.seqera.io/library/portello:0.7.0--e30f230f4d2812dd'}"

    input:
    tuple val(meta), path(asm_to_ref_bam), path(asm_to_ref_bai), path(read_to_asm_bam), path(read_to_asm_bai), path(ref_fasta), val(assembly_mode), val(output_vcf)

    output:
    tuple val(meta), path("*_remapped.bam"), emit: bam
    tuple val(meta), path("*_unassembled.bam"), emit: unassembled_bam
    tuple val(meta), path("*.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val("${task.process}"), val('portello'), eval("portello --version | sed -e 's/portello //'"), topic: versions, emit: versions_portello

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_output = output_vcf ? "--phased-het-vcf-prefix ${prefix}" : ""
    """
    portello \
        ${args} \
        --threads ${task.cpus} \
        --ref ${ref_fasta} \
        --assembly-to-ref ${asm_to_ref_bam} \
        --read-to-assembly ${read_to_asm_bam} \
        --input-assembly-mode ${assembly_mode} \
        --unassembled-read-output ${prefix}_unassembled.bam \
        ${vcf_output} \
        --remapped-read-output ${prefix}_remapped.bam
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_output = output_vcf ? "echo | gzip -c > ${prefix}.vcf.gz; touch ${prefix}.vcf.gz.tbi" : ''
    """
    echo ${args}

    ${vcf_output}
    touch ${prefix}_unassembled.bam
    touch ${prefix}_remapped.bam
    """
}

process BISCUIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0':
        'biocontainers/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0' }"

    input:
    tuple val(meta), path(normal_bams), path(normal_bais), path(tumor_bam), path(tumor_bai)
    path index

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biscuit_cpus = (int) Math.max(Math.floor(task.cpus*0.9),1)
    def bgzip_cpus = task.cpus-biscuit_cpus
    if ( tumor_bam != [] && normal_bams.toList().size() > 1 ) error "[BISCUIT_PILEUP] error: Tumor BAM provided with more than one normal BAM"
    if ( tumor_bam.toList().size() > 1 ) error "[BISCUIT_PILEUP] error: more than one tumor BAM provided"
    input = ( tumor_bam==[] ) ? "${normal_bams}" : "-S -T ${tumor_bam} -I ${normal_bams}"
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/\\.bis.amb\$//'`

    biscuit pileup \\
        -@ $biscuit_cpus \\
        $args \\
        \$INDEX \\
        $input \\
        | bgzip -@ $bgzip_cpus $args2 > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """
}

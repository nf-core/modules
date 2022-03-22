process BISCUIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113 bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0':
        'quay.io/biocontainers/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0' }"

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
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

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
}

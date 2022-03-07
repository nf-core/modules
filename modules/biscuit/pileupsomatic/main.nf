process BISCUIT_PILEUPSOMATIC {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113 bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0':
        'quay.io/biocontainers/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0' }"

    input:
    tuple val(meta), path(normal_bam), path(tumor_bam), path(normal_bai), path(tumor_bai)
    path index

    output:
    tuple val(meta), path("*.somatic.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biscuit_cpus = (int) Math.max(Math.floor(task.cpus*0.9),1)
    def bgzip_cpus = task.cpus-biscuit_cpus
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

    biscuit pileup \\
        -@ $biscuit_cpus \\
        $args \\
        \$INDEX \\
        -S \\
        -I $normal_bam \\
        -T $tumor_bam \\
        | bgzip -@ $bgzip_cpus $args2 > ${prefix}.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """
}

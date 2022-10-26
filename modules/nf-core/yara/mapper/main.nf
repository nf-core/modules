process YARA_MAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::yara=1.0.2 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f13549097a0d1ca36f9d4f017636fb3609f6c083:d6c969c1e20cc02a9234961c07a24bb0887f05ea-0' :
        'quay.io/biocontainers/mulled-v2-f13549097a0d1ca36f9d4f017636fb3609f6c083:d6c969c1e20cc02a9234961c07a24bb0887f05ea-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.mapped.bam"), emit: bam
    tuple val(meta), path("*.mapped.bam.bai"), emit: bai
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_prefix = index[0].baseName.substring(0,index[0].baseName.lastIndexOf('.'))
    if (meta.single_end) {
        """
        yara_mapper \\
            $args \\
            -t $task.cpus \\
            -f bam \\
            ${index_prefix} \\
            $reads | samtools view -@ $task.cpus -hb -F4 | samtools sort -@ $task.cpus > ${prefix}.mapped.bam

        samtools index -@ $task.cpus ${prefix}.mapped.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            yara: \$(echo \$(yara_mapper --version 2>&1) | sed 's/^.*yara_mapper version: //; s/ .*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        yara_mapper \\
            $args \\
            -t $task.cpus \\
            -f bam \\
            ${index_prefix} \\
            ${reads[0]} \\
            ${reads[1]} > output.bam

        samtools view -@ $task.cpus -hF 4 -f 0x40 -b output.bam | samtools sort -@ $task.cpus > ${prefix}_1.mapped.bam
        samtools view -@ $task.cpus -hF 4 -f 0x80 -b output.bam | samtools sort -@ $task.cpus > ${prefix}_2.mapped.bam

        samtools index -@ $task.cpus ${prefix}_1.mapped.bam
        samtools index -@ $task.cpus ${prefix}_2.mapped.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            yara: \$(echo \$(yara_mapper --version 2>&1) | sed 's/^.*yara_mapper version: //; s/ .*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}

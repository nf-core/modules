process YARA_MAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f13549097a0d1ca36f9d4f017636fb3609f6c083:de7982183b85634270540ac760c2644f16e0b6d1-0' :
        'quay.io/biocontainers/mulled-v2-f13549097a0d1ca36f9d4f017636fb3609f6c083:de7982183b85634270540ac760c2644f16e0b6d1-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.mapped.bam")    , emit: bam
    tuple val(meta), path("*.mapped.bam.bai"), emit: bai
    tuple val("${task.process}"), val('yara'), eval("yara_mapper --version 2>&1 | grep 'yara_mapper version' | sed 's/^.*yara_mapper version: //; s/ .*\$//'"), topic: versions, emit: versions_yara
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | head -n1 | sed 's/^.*samtools //'"), topic: versions, emit: versions_samtools


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
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}.mapped.bam
        touch ${prefix}.mapped.bam.bai
        """
    } else {
        """
        touch ${prefix}_1.mapped.bam
        touch ${prefix}_2.mapped.bam.bai
        """
    }

}

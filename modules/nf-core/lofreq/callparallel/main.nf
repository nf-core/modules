process LOFREQ_CALLPARALLEL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4' :
        'biocontainers/lofreq:2.1.5--py38h588ecb2_4' }"

    input:
    tuple val(meta), path(tumor), path(tumor_index), path(intervals)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_intervals = intervals ? "-l ${intervals}" : ""

    def tumor_cram =  tumor.Extension == "cram" ? true : false
    def tumor_bam = tumor.Extension == "bam" ? true : false
    def tumor_out = tumor_cram ? tumor.BaseName + ".bam" : "${tumor}"

    def samtools_cram_convert = ''
    samtools_cram_convert += tumor_cram ? "    samtools view -T $fasta $tumor -@ $task.cpus -o $tumor_out\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools index $tumor_out\n" : ''
    """
    $samtools_cram_convert

    lofreq \\
        call-parallel \\
        --pp-threads $task.cpus \\
        $args \\
        $options_intervals \\
        -f $fasta \\
        -o ${prefix}.vcf.gz \\
        $tumor_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}

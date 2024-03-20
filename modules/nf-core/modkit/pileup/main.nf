process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::ont-modkit=0.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.2.2--hdcf5f25_0':
        'biocontainers/ont-modkit:0.2.2--hdcf5f25_0' }"

    input:
    tuple val(meta), path(bam)
    path bai
    path out_bed

    output:
    tuple val(meta), path(out_bed), emit: bed
    tuple val(meta), path("*_pileup.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filename = "${bam}"[0..<"${bam}".lastIndexOf('.')]
    """
    mkdir -p tmp/ &&
    cp $bam $bai tmp/ &&
    modkit \\
        pileup \\
        tmp/$bam \\
        $out_bed \\
        --threads $task.cpus \\
        --log-filepath ${filename}_pileup.log \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(modkit --version 2>&1) | sed 's/^.*modkit //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filename = "${bam}"[0..<"${bam}".lastIndexOf('.')]
    """
    touch $out_bed
    touch ${filename}_pileup.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

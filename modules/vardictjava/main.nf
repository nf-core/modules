process VARDICTJAVA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vardict-java=1.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0':
        'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(regions_of_interest)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vardict-java \\
        $args \\
        -c 1 -S 2 -E 3 \\
        -b $bam \\
        -th $task.cpus \\
        -N $prefix \\
        -G $reference_fasta \\
        $regions_of_interest \\
        | gzip -c > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: \$(echo 1.8.3)
    END_VERSIONS
    """
}

process GRAPHTYPER_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::graphtyper=2.7.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphtyper:2.7.2--h7d7f7ad_0':
        'quay.io/biocontainers/graphtyper:2.7.2--h7d7f7ad_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path ref_fai
    path region

    output:
    tuple val(meta), path("results/*/*.vcf.gz"), emit: vcf
    tuple val(meta), path("results/*/*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_path_text = bam.join('\n')
    """
    printf "$bam_path_text" > bam_list.txt
    graphtyper \\
        genotype \\
        $args \\
        $ref \\
        --threads $task.cpus \\
        --sams bam_list.txt \\
        --region_file $region

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """
}

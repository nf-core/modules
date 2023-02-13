process GRAPHTYPER_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::graphtyper=2.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphtyper:2.7.2--h7d7f7ad_0':
        'quay.io/biocontainers/graphtyper:2.7.2--h7d7f7ad_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(ref_fai)
    path region_file  // can be empty if --region is supplied to task.ext.args

    output:
    tuple val(meta), path("results/*/*.vcf.gz"), emit: vcf
    tuple val(meta), path("results/*/*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_path_text = bam.join('\\n')
    def region_text = region_file.size() > 0 ? "--region_file ${region_file}" : ""
    if (region_file.size() == 0 && ! args.contains("region")) {
        error "GRAPHTYPER_GENOTYPE requires either a region file or a region specified using '--region' in ext.args"
    }
    """
    printf "$bam_path_text" > bam_list.txt
    graphtyper \\
        genotype \\
        $ref \\
        $args \\
        --sams bam_list.txt \\
        --threads $task.cpus \\
        $region_text

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """
}

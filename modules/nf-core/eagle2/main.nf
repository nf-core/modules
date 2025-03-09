process EAGLE2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eagle2:2.4.1--h6a68c12_0':
        'biocontainers/eagle2:2.4.1--h6a68c12_0' }"

    input:
    tuple val(meta), path(input), path(index), path(ref_vcf), path(ref_index), path(map)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf}"), emit: phased_variants
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_cmd = ref_vcf ? "--vcfTarget $input" : "--vcf $input"
    def vcfref_cmd = ref_vcf ? "--vcfRef $ref_vcf" : ""
    def remove_cmd = remove ? "--remove $remove" : ""
    def exclude_cmd = exclude ? "--exclude $exclude" : ""
    """
    eagle \\
        $input_cmd \\
        $vcfref_cmd \\
        --geneticMapFile $map \\
        $remove_cmd \\
        $exclude_cmd \\
        $args \\
        --numThreads $task.cpus \\
        --outPrefix ${prefix} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eagle2: \$(eagle --help | sed -n 's/.*Eagle v\\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--vcfOutFormat z") ? "vcf.gz" :
                    args.contains("--vcfOutFormat v") ? "vcf"    :
                    args.contains("--vcfOutFormat b") ? "bcf"    :
                    args.contains("--vcfOutFormat u") ? "bcf"    :
                    "vcf.gz"
    def create_cmd = extension.endsWith("gz") ? "echo '' | bgzip >" : "touch"
    """
    ${create_cmd} ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eagle2: \$(eagle --help | sed -n 's/.*Eagle v\\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}


process SVTYPER_SVTYPER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::svtyper=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtyper:0.7.0--py27h24bf2e0_1':
        'quay.io/biocontainers/svtyper:0.7.0--py27h24bf2e0_1' }"

    input:
    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.json")  , emit: bam
    tuple val(meta), path("*.gt.vcf"), emit: gt_vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf  = vcf ? "--input_vcf ${vcf}" : ""

    """
    svtyper \\
        $vcf \\
        --bam $bam \\
        --lib_info ${prefix}.json \\
        --output_vcf ${prefix}.gt.vcf \\
        --ref_fasta $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtyper: \$(echo \$(svtyper) | sed 's/[^0-9.]*\\([0-9.]*\\).*/\\1/' )
    END_VERSIONS
    """
}

process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.0.10--py39hd649744_1':
        'biocontainers/clair3:1.0.10--py39hd649744_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(model),val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*merge_output.vcf.gz"),            emit: vcf
    tuple val(meta), path("*merge_output.vcf.gz.tbi"),        emit: tbi
    tuple val(meta), path("*phased_merge_output.vcf.gz"),     emit: phased_vcf, optional: true
    tuple val(meta), path("*phased_merge_output.vcf.gz.tbi"), emit: phased_tbi, optional: true
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clair3.sh \\
        --bam_fn=$bam \\
        --ref_fn=$reference \\
        --threads=$task.cpus \\
        --output=\$PWD \\
        --platform=$platform \\
        --model=$model \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh  --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.phased_merge_output.vcf.gz
    echo "" | gzip > ${prefix}.phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}.merge_output.vcf.gz
    echo "" | gzip > ${prefix}.merge_output.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """
}

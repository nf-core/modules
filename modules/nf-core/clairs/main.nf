process CLAIRS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/hkubal/clairs:v0.4.1':
        'docker.io/hkubal/clairs:v0.4.1' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(model)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.vcf.gz"),      emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"),  emit: tbi
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /opt/bin/run_clairs \
        --tumor_bam_fn $tumor_bam \\
        --normal_bam_fn $normal_bam \\
        --ref_fn $reference \\
        --threads $task.cpus \\
        --platform $model \\
        --output_dir . \\
        --output_prefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(/opt/bin/run_clairs --version |& sed '1!d; s|run_clairs ||'
)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(/opt/bin/run_clairs --version |& sed '1!d; s|run_clairs ||'
)
    END_VERSIONS
    """
}

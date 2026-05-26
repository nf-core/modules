
process CLAIRSTO {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/hkubal/clairs-to:v0.4.0':
        'docker.io/hkubal/clairs-to:v0.4.0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), val(model)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*/indel.vcf.gz"),     emit: indel_vcf
    tuple val(meta), path("*/indel.vcf.gz.tbi"),     emit: indel_tbi
    tuple val(meta), path("*/snv.vcf.gz"),       emit: snv_vcf
    tuple val(meta), path("*/snv.vcf.gz.tbi"),       emit: snv_tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def conda_prefix = workflow.containerEngine == 'singularity' ? '--conda_prefix /opt/micromamba/envs/clairs-to' : ''

    """
    /opt/bin/run_clairs_to \
        --tumor_bam_fn $tumor_bam \\
        --ref_fn $reference \\
        --platform $model \\
        --threads $task.cpus \\
        --output_dir $prefix \\
        $conda_prefix \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairsto: \$(/opt/bin/run_clairs_to --version |& sed '1!d ; s/run_clairs_to //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p output
    echo "" | gzip > output/snv.vcf.gz
    touch output/snv.vcf.gz.tbi
    echo "" | gzip > output/indel.vcf.gz
    touch output/indel.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairsto: \$(/opt/bin/run_clairs_to --version |& sed '1!d ; s/run_clairs_to //')
    END_VERSIONS
    """
}

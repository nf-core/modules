
process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clair3:v1.0.10':
        'hkubal/clair3:v1.0.10' }"

    input:
    tuple val(meta), path(bam), path(bai), path(model), val(platform)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*merge_output.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def clair3_command = workflow.containerEngine ? "/opt/bin/run_clair3.sh" : "run_clair3.sh"

    """
    $clair3_command \\
        --bam_fn=$bam \\
        --ref_fn=$fasta \\
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
    touch ${prefix}.phased_merge_output.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(/opt/bin/run_clair3.sh --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """
}

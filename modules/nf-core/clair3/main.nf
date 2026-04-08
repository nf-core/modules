process CLAIR3 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://hkubal/clair3:v2.0.0' :
    'docker.io/hkubal/clair3:v2.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(packaged_model), path(user_model), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("${prefix}merge_output.vcf.gz"),            emit: vcf
    tuple val(meta), path("${prefix}merge_output.vcf.gz.tbi"),        emit: tbi
    tuple val(meta), path("${prefix}phased_merge_output.vcf.gz"),     emit: phased_vcf, optional: true
    tuple val(meta), path("${prefix}phased_merge_output.vcf.gz.tbi"), emit: phased_tbi, optional: true
    tuple val(meta), path("${prefix}merge_output.gvcf.gz"),           emit: gvcf, optional: true
    tuple val(meta), path("${prefix}merge_output.gvcf.gz.tbi"),       emit: gtbi, optional: true
    tuple val("${task.process}"), val('clair3'), eval('run_clair3.sh --version | sed "s/^Clair3 v//"'), emit: versions_clair3, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def model = ""
    if (!user_model) {
        if (workflow.containerEngine in ['singularity', 'docker', 'podman']) {
            model = "/opt/models/${packaged_model}"
        } else {
            error "Clair3 packaged models are only available in Docker/Singularity/Podman containers. Please use one of these profiles or provide a user_model instead."
        }
    }
    if (!packaged_model) {
        model = "${user_model}"
    }
    if (packaged_model && user_model) {
        error("Two models specified ${user_model} and ${packaged_model}, specify one of them.")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clair3.sh \\
        --bam_fn=\${PWD}/${bam} \\
        --ref_fn=\${PWD}/${reference} \\
        --threads=${task.cpus} \\
        --output=. \\
        --platform=${platform} \\
        --model_path=${model} \\
        ${args}

    # Rename to add prefix
    for file in merge_output.vcf.gz \
            merge_output.vcf.gz.tbi \
            phased_merge_output.vcf.gz \
            phased_merge_output.vcf.gz.tbi \
            merge_output.gvcf.gz \
            merge_output.gvcf.gz.tbi; do
    if [ -e "\$file" ]; then
        mv "\$file" "${prefix}\$file"
    fi
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}phased_merge_output.vcf.gz
    touch ${prefix}phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}merge_output.vcf.gz
    touch ${prefix}merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}merge_output.gvcf.gz
    touch ${prefix}merge_output.gvcf.gz.tbi
    """
}

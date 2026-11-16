// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CLAIR3 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fac3c639be51f2aa2d33b51009f2a6c4175392b90c1073c670f7e70afa367f6d/data'
        : 'community.wave.seqera.io/library/clair3:2.0.5--af64dfa2c8f0b388'}"

    input:
    tuple val(meta), path(bam), path(bai), val(packaged_model), path(user_model), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)
    path(bed)

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
    if (packaged_model && user_model) {
        error("Two models specified ${user_model} and ${packaged_model}, specify one of them.")
    }
    // bioconda clair3 bundles its models next to run_clair3.sh, wherever the environment is installed
    def model = user_model ? "${user_model}" : "\$(dirname \$(command -v run_clair3.sh))/models/${packaged_model}"
    def args = task.ext.args ?: ''
    def bed_arg = bed ? "--bed_fn=\${PWD}/${bed}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clair3.sh \\
        --bam_fn=\${PWD}/${bam} \\
        --ref_fn=\${PWD}/${reference} \\
        --threads=${task.cpus} \\
        --output=. \\
        --platform=${platform} \\
        --model_path=${model} \\
        ${bed_arg} \\
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

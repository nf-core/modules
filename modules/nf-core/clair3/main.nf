process CLAIR3 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e70b0f4389028f4dc88efde1aac7139927c898cf7add680e14724d97fecd3d32/data'
        : 'community.wave.seqera.io/library/clair3:1.2.0--b1b03d4e9d1b6a2e'}"

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
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def model = ""
    if (!user_model) {
        // In seqera containers `MAMBA_ROOT_PREFIX` is always available
        // whereas `CONDA_PREFIX` may not be
        // see https://github.com/seqeralabs/wave/issues/886
        model = "\${CONDA_PREFIX:-\$MAMBA_ROOT_PREFIX}/bin/models/${packaged_model}"
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
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --output=. \\
        --platform=${platform} \\
        --model=${model} \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh  --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}phased_merge_output.vcf.gz
    touch ${prefix}phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}merge_output.vcf.gz
    touch ${prefix}merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}merge_output.gvcf.gz
    touch ${prefix}merge_output.gvcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """
}

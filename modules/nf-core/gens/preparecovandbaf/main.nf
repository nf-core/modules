process PREPAREGENSINPUTDATA {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/2931cf8e5eac90d974b42ec069b4268fd7f8aefa37b46fa18eb342a4802a967d/data':
        'community.wave.seqera.io/library/tabix_pip_gens-input-data-tools:acc3fd1b79233d2d' }"

    input:
    tuple val(meta), path(read_counts), path(gvcf), path(gvcf_tbi)
    path baf_positions

    output:
    tuple val(meta), path("*.cov.bed.gz")     , emit: cov_gz
    tuple val(meta), path("*.cov.bed.gz.tbi") , emit: cov_tbi
    tuple val(meta), path("*.baf.bed.gz")     , emit: baf_gz
    tuple val(meta), path("*.baf.bed.gz.tbi") , emit: baf_tbi
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "The gens pre-processing module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // FIXME: Can I avoid this
    def python_base = "/opt/conda/lib/python3.14/site-packages/gens_input_data_tools"
    """
    python3 $python_base/generate_cov_and_baf.py \\
        --coverage read_counts \\
        --gvcf $gvcf \\
        --label $prefix \\
        --baf_positions $baf_positions \\
        --bgzip_tabix_output \\
        $args \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preparegensinputdata: \$(python3 $python_base/generate_cov_and_baf.py --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def python_base = "/opt/conda/lib/python3.14/site-packages/gens_input_data_tools"
    """
    echo "" | gzip > ${prefix}.cov.bed.gz
    touch ${prefix}.cov.bed.gz.tbi
    echo "" | gzip > ${prefix}.baf.bed.gz
    touch ${prefix}.baf.bed.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preparegensinputdata: \$(python3 $python_base/generate_cov_and_baf.py --version)
    END_VERSIONS
    """
}

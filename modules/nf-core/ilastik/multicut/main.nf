process ILASTIK_MULTICUT {
    tag "$meta.id"
    label 'process_low'


    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "ILASTIK_MULTICUT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "biocontainers/ilastik:1.4.0_cv1"

    input:
    tuple val(meta), path(h5)
    tuple val(meta2), path (ilp)
    tuple val(meta3), path (probs)

    output:
    tuple val(meta), path("*.tiff") , emit: out_tiff
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /opt/ilastik-1.4.0-Linux/run_ilastik.sh \\
        --headless \\
        --readonly 1 \\
        --project=$ilp \\
        --raw_data=$h5 \\
        --probabilities=$probs \\
        --export_source="Multicut Segmentation" \\
        --output_filename_format=${prefix}.tiff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(/opt/ilastik-1.4.0-Linux/run_ilastik.sh --headless --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.4.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.tiff
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik:: $VERSION
    END_VERSIONS
    """
}

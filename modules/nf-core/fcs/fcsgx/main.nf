process FCS_FCSGX {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-gx.sif':
        'docker.io/ncbi/fcs-gx:0.4.0' }"

    input:
    tuple val(meta), path(assembly)
    path gxdb

    output:
    tuple val(meta), path("out/*.fcs_gx_report.txt"), emit: fcs_gx_report
    tuple val(meta), path("out/*.taxonomy.rpt")     , emit: taxonomy_report
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/fcsgx/rungx

Reason:
This module is now renamed as FCSGX_RUNGX and as been updated to the latest version
"""
    // Comment out this block to disable the deprecation warning
    assert false: deprecation_message

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def FCSGX_VERSION = '0.4.0'

    """
    python3 /app/bin/run_gx \\
        --fasta $assembly \\
        --out-dir ./out \\
        --gx-db ./${gxdb[0].baseName} \\
        --tax-id ${meta.taxid} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """

    stub:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/fcsgx/rungx

Reason:
This module is now renamed as FCSGX_RUNGX and as been updated to the latest version
"""
    // Comment out this block to disable the deprecation warning
    assert false: deprecation_message

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.4.0'

    """
    mkdir -p out
    touch out/${prefix}.fcs_gx_report.txt
    touch out/${prefix}.taxonomy.rpt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """
}

process FCS_FCSGX {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.2.3/fcs-gx.0.2.3.sif':
        'ncbi/fcs-gx:0.2.3' }"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(assembly)
    path gxdb

    output:
    tuple val(meta), path("out/*.fcs_gx_report.txt"), emit: fcs_gx_report
    tuple val(meta), path("out/*.taxonomy.rpt")     , emit: taxonomy_report
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.2.3'

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
}

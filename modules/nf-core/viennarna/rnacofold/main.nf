process VIENNARNA_RNACOFOLD {
    tag "$rnacofold_fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9fa2645e08f7c6ed723b1548b6782740bf59fe79d62046bb8fb0ec6326d2525f/data'
:         'community.wave.seqera.io/library/viennarna:2.6.4--bd1078d77d380b29' }"

    input:
    tuple val(meta), path(rnacofold_fasta)

    output:
    tuple val(meta), path("*.csv")                                                                   , emit: rnacofold_csv
    tuple val(meta), path("*.ps")                                                                    , emit: rnacofold_ps
    tuple val("${task.process}"), val('RNAcofold'), eval("RNAcofold --version 2>&1 | sed -n '1s/RNAcofold //p'"), emit: versions_rnacofold, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNAcofold < ${rnacofold_fasta} \\
        --jobs=${task.cpus} \\
        --output-format="D" \\
        ${args} \\
        > ${rnacofold_fasta.baseName}.csv
    """

    stub:
    """
    touch ${rnacofold_fasta.baseName}.csv
    touch ${rnacofold_fasta.baseName}.ps
    """
}

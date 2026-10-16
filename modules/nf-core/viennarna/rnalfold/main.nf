process VIENNARNA_RNALFOLD {
    tag '$rna_fastq'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9fa2645e08f7c6ed723b1548b6782740bf59fe79d62046bb8fb0ec6326d2525f/data'
:         'community.wave.seqera.io/library/viennarna:2.6.4--bd1078d77d380b29' }"

    input:
    path fasta

    output:
    path "*.lfold"                                                                       , emit: rnalfold_txt
    tuple val("${task.process}"), val('RNALfold'), eval("RNALfold --version 2>&1 | sed -n '1s/RNALfold //p'"), emit: versions_rnalfold, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNALfold \\
        ${args} \\
        --infile=$fasta \\
        --outfile=${fasta.baseName}.lfold
    """

    stub:
    """
    touch ${fasta.baseName}.lfold
    touch ${fasta}.ps
    """
}

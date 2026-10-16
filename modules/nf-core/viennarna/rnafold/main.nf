process VIENNARNA_RNAFOLD {
    tag '$rna_fastq'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9fa2645e08f7c6ed723b1548b6782740bf59fe79d62046bb8fb0ec6326d2525f/data'
:         'community.wave.seqera.io/library/viennarna:2.6.4--bd1078d77d380b29' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fold")                                                             , emit: rnafold_txt
    tuple val(meta), path("*.ps")                                                               , emit: rnafold_ps
    tuple val("${task.process}"), val('RNAfold'), eval("RNAfold --version 2>&1 | sed -n '1s/RNAfold //p'"), emit: versions_rnafold, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id ?: "${fasta.getName()}"
    """
    RNAfold \\
        ${args} \\
        --jobs=${task.cpus} \\
        --infile=${fasta} \\
        --outfile=${prefix}.fold
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id ?: "${fasta.getName()}"
    """
    touch ${prefix}.fold
    touch ${prefix}.ps
    """
}

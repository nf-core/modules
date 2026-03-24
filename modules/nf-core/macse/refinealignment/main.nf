process MACSE_REFINEALIGNMENT {
    tag "${meta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/macse:2.07--hdfd78af_0'
        : 'biocontainers/macse:2.07--hdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_NT.{fa,fas,fasta,aln}"), emit: nt
    tuple val(meta), path("*_AA.{fa,fas,fasta,aln}"), emit: aa
    tuple val("${task.process}"), val("macse"), eval("macse --version | sed -n 's/.*V\\([0-9]*\\.[0-9]*\\).*/\\1/p'"), topic: versions, emit: versions_macse

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName

    """
    sed '/^>/!s/[[:space:]]//g' ${fasta} > ${prefix}.clean.fasta

    macse -prog refineAlignment \\
        -align ${prefix}.clean.fasta \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    echo ${args}
    touch ${prefix}_NT.fas
    touch ${prefix}_AA.fas
    """
}

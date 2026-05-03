process RAPIDNJ {
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-805c6e0f138f952f9c61cdd57c632a1a263ea990:3c52e4c8da6b3e4d69b9ca83fa4d366168898179-0' :
        'quay.io/biocontainers/mulled-v2-805c6e0f138f952f9c61cdd57c632a1a263ea990:3c52e4c8da6b3e4d69b9ca83fa4d366168898179-0' }"

    input:
    path alignment

    output:
    path "*.sth"       , emit: stockholm_alignment
    path "*.tre"       , emit: phylogeny
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('rapidnj'), eval('echo 2.3.2'), emit: versions_rapidnj, topic: versions
    tuple val("${task.process}"), val('biopython'), eval('python -c "import Bio; print(Bio.__version__)"'), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python \\
        -c 'from Bio import SeqIO; SeqIO.convert("$alignment", "fasta", "alignment.sth", "stockholm")'

    rapidnj \\
        alignment.sth \\
        $args \\
        -i sth \\
        -c $task.cpus \\
        -x rapidnj_phylogeny.tre
    """

    stub:
    """
    touch alignment.sth
    touch rapidnj_phylogeny.tre
    """
}

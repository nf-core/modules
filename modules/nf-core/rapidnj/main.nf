process RAPIDNJ {
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::rapidnj=2.3.2 conda-forge::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-805c6e0f138f952f9c61cdd57c632a1a263ea990:3c52e4c8da6b3e4d69b9ca83fa4d366168898179-0' :
        'biocontainers/mulled-v2-805c6e0f138f952f9c61cdd57c632a1a263ea990:3c52e4c8da6b3e4d69b9ca83fa4d366168898179-0' }"

    input:
    path alignment

    output:
    path "*.sth"       , emit: stockholm_alignment
    path "*.tre"       , emit: phylogeny
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '2.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    python \\
        -c 'from Bio import SeqIO; SeqIO.convert("$alignment", "fasta", "alignment.sth", "stockholm")'

    rapidnj \\
        alignment.sth \\
        $args \\
        -i sth \\
        -c $task.cpus \\
        -x rapidnj_phylogeny.tre

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rapidnj: $VERSION
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

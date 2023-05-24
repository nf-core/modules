process GENOMAD_ENDTOEND {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::genomad=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.5.2--pyhdfd78af_0':
        'biocontainers/genomad:1.5.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path genomad_db

    output:
    tuple val(meta), path("*_summary/*_plasmid.fna")            , emit: plasmid_fasta
    tuple val(meta), path("*_summary/*_plasmid_genes.tsv")      , emit: plasmid_genes
    tuple val(meta), path("*_summary/*_plasmid_proteins.faa")   , emit: plasmid_proteins
    tuple val(meta), path("*_summary/*_plasmid_summary.tsv")    , emit: plasmid_summary
    tuple val(meta), path("*_summary/*_virus.fna")              , emit: virus_fasta
    tuple val(meta), path("*_summary/*_virus_genes.tsv")        , emit: virus_genes
    tuple val(meta), path("*_summary/*_virus_proteins.faa")     , emit: virus_proteins
    tuple val(meta), path("*_summary/*_virus_summary.tsv")      , emit: virus_summary
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    genomad \\
        end-to-end \\
        $fasta \\
        ./ \\
        $genomad_db \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}

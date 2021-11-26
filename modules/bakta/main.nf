process BAKTA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bakta=1.2.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.2.2--pyhdfd78af_0' :
        'quay.io/biocontainers/bakta:1.2.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}.embl")             , emit: embl
    tuple val(meta), path("${prefix}.faa")              , emit: faa
    tuple val(meta), path("${prefix}.ffn")              , emit: ffn
    tuple val(meta), path("${prefix}.fna")              , emit: fna
    tuple val(meta), path("${prefix}.gbff")             , emit: gbff
    tuple val(meta), path("${prefix}.gff3")             , emit: gff
    tuple val(meta), path("${prefix}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}.tsv")              , emit: tsv
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    """
    bakta \\
        $args \\
        --threads $task.cpus \\
        --prefix $prefix \\
        --db $db \\
        $proteins_opt \\
        $prodigal_tf \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    touch ${prefix}.embl
    touch ${prefix}.faa
    touch ${prefix}.ffn
    touch ${prefix}.fna
    touch ${prefix}.gbff
    touch ${prefix}.gff3
    touch ${prefix}.hypotheticals.tsv
    touch ${prefix}.hypotheticals.faa
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}

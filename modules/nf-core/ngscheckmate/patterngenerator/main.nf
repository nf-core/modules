process NGSCHECKMATE_PATTERNGENERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::ngscheckmate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' :
        'biocontainers/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(fasta)
    path(bowtie_index)

    output:
    tuple val(meta), path("*.pt"), emit: pt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "makesnvpattern.pl generates a fasta file with the same name as the input fasta, use \"task.ext.prefix\" to disambiguate!"

    """
    makesnvpattern.pl ${bed} ${fasta} ${bowtie_index}/${fasta.getBaseName()} . ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "makesnvpattern.pl generates a fasta file with the same name as the input fasta, use \"task.ext.prefix\" to disambiguate!"

    """
    touch ${prefix}.pt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """
}

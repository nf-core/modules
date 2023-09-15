process NGSCHECKMATE_PATTERNGENERATOR {
    tag "$meta1.id"
    label 'process_single'

    conda "bioconda::ngscheckmate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' :
        'biocontainers/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' }"

    input:
    tuple val(meta1), path(bed)
    tuple val(meta2), path(fasta)
    path(bowtie_index)

    output:
    tuple val(meta1), path("*.pt"), emit: pt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta1.id}"

    """

    mkdir -p outdir
    makesnvpattern.pl ${bed} ${fasta} ${bowtie_index}/${fasta.getBaseName()} outdir ${prefix}
    cp outdir/${prefix}.pt .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """


    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """
}

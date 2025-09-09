process TRANSDECODER_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.7.1--pl5321hdfd78af_0':
        'biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(fold)

    output:
    tuple val(meta), path("*.transdecoder.pep")  , emit: pep
    tuple val(meta), path("*.transdecoder.gff3") , emit: gff3
    tuple val(meta), path("*.transdecoder.cds")  , emit: cds
    tuple val(meta), path("*.transdecoder.bed")  , emit: bed
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    TransDecoder.Predict \\
        $args \\
        -O . \\
        -t \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.Predict --version) | sed -e "s/TransDecoder.Predict //g")
    END_VERSIONS
    """

    stub:
    def fasta_no_gz = fasta.toString() - '.gz'
    """
    touch ${fasta_no_gz}.transdecoder.pep
    touch ${fasta_no_gz}.transdecoder.gff3
    touch ${fasta_no_gz}.transdecoder.cds
    touch ${fasta_no_gz}.transdecoder.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.Predict --version) | sed -e "s/TransDecoder.Predict //g")
    END_VERSIONS
    """
}

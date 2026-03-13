process TRANSDECODER_LONGORF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/transdecoder:5.7.1--pl5321hdfd78af_0' :
    'biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${output_dir_name}/*.pep")   , emit: pep
    tuple val(meta), path("${output_dir_name}/*.gff3")  , emit: gff3
    tuple val(meta), path("${output_dir_name}/*.cds")   , emit: cds
    tuple val(meta), path("${output_dir_name}/*.dat")   , emit: dat
    path("${output_dir_name}")                          , emit: folder
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args         ?: ''
    def prefix      = task.ext.prefix       ?: "${meta.id}"
    def fasta_no_gz = fasta.toString()      - '.gz'
    output_dir_name = "${meta.id}/${fasta_no_gz}.transdecoder_dir"
    """
    TransDecoder.LongOrfs \\
        $args \\
        -O $prefix \\
        -t \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.LongOrfs --version) | sed -e "s/TransDecoder.LongOrfs //g")
    END_VERSIONS
    """

    stub:
    def fasta_no_gz = fasta.toString()      - '.gz'
    output_dir_name = "${meta.id}/${fasta_no_gz}.transdecoder_dir"
    """
    mkdir -p $output_dir_name

    touch ${output_dir_name}/longest_orfs.pep
    touch ${output_dir_name}/longest_orfs.gff3
    touch ${output_dir_name}/longest_orfs.cds
    touch ${output_dir_name}/base_freqs.dat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.LongOrfs --version) | sed -e "s/TransDecoder.LongOrfs //g")
    END_VERSIONS
    """
}

process DASTOOL_SCAFFOLDS2BIN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::das_tool=1.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/das_tool:1.1.3--r41hdfd78af_0' :
        'quay.io/biocontainers/das_tool:1.1.3--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    val(extension)

    output:
    tuple val(meta), path("*.tsv"), emit: scaffolds2bin
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = extension ? extension : "fasta"
    def clean_fasta = fasta.toString() - ".gz"
    def decompress_fasta = fasta.toString() == clean_fasta ? "" : "gunzip -q -f $fasta"
    """
    $decompress_fasta

    Fasta_to_Scaffolds2Bin.sh \\
        $args \\
        -i . \\
        -e $file_extension \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dastool: \$( DAS_Tool --version 2>&1 | grep "DAS Tool" | sed 's/DAS Tool version //' )
    END_VERSIONS
    """
}

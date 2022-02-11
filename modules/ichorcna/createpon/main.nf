def VERSION = '0.3.2' // Version information not provided by tool on CLI

process ICHORCNA_CREATEPON {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::r-ichorcna=0.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.3.2--r41hdfd78af_0' :
        'quay.io/biocontainers/r-ichorcna:0.3.2--r41hdfd78af_0' }"

    input:
    path wigs
    path gc_wig
    path map_wig
    path centromere

    output:
    path "*.rds"        , emit: rds
    path "*.txt"        , emit: txt
    path "versions.yml" , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def prefix = task.ext.prefix ?: "PoN"

    """
    echo ${wigs} | tr " " "\\n" > wig_files.txt

    createPanelOfNormals.R \\
        --filelist wig_files.txt \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${centro} \\
        ${args} \\
        --outfile ${prefix}

    rm wig_files.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
}

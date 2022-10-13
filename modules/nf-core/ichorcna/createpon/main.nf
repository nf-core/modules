process ICHORCNA_CREATEPON {
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::r-ichorcna=0.3.2" : null)
    def container_image = "/r-ichorcna:0.3.2--r41hdfd78af_0"
                                                     container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

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
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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

process ICHORCNA_RUN {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.5.0--pl5321r42hdfd78af_0' :
        'docker.io/scwatts/r-ichorcna:0.5.1--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(wig)
    path gc_wig
    path map_wig
    path normal_wig
    path normal_background
    path centromere
    path rep_time_wig
    path exons

    output:
    tuple val(meta), path("${prefix}.RData")             , emit: rdata
    tuple val(meta), path("${prefix}.seg")               , emit: seg
    tuple val(meta), path("${prefix}.cna.seg")           , emit: cna_seg
    tuple val(meta), path("${prefix}.seg.txt")           , emit: seg_txt
    tuple val(meta), path("${prefix}.correctedDepth.txt"), emit: corrected_depth
    tuple val(meta), path("${prefix}.params.txt")        , emit: ichorcna_params
    tuple val(meta), path("${prefix}/*.pdf")             , emit: plots
    tuple val(meta), path("**/${prefix}_genomeWide.pdf") , emit: genome_plot
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args       ?: ''
    prefix = task.ext.prefix   ?: "${meta.id}"
    def norm   = normal_wig        ? "normal_wig='${normal_wig}',"          : 'normal_wig=NULL,'
    def pon    = normal_background ? "normal_panel='${normal_background}'," : 'normal_panel=NULL,'
    def map    = map_wig           ? "mapWig='${map_wig}',"                 : 'mapWig=NULL,'
    def centro = centromere        ? "centromere='${centromere}',"          : ''
    def rep    = rep_time_wig      ? "repTimeWig='${rep_time_wig}',"        : 'repTimeWig=NULL,'
    def exons  = exons             ? "exons.bed='${exons}',"                : ''
    """
    #!/usr/bin/env Rscript
    library("ichorCNA")
    library("yaml")

    run_ichorCNA(
        tumor_wig='${wig}',
        id='${prefix}',
        cores=${task.cpus},
        gcWig='${gc_wig}',
        $norm
        $pon
        $map
        $centro
        $rep
        $exons
        $args
        outDir="."
    )


    ### Make Versions YAML for NF-Core ###
    versions = list()
    versions["r"]        <- paste(R.Version()\$major, R.Version()\$minor, sep=".")
    versions["ichorCNA"] <- paste(packageVersion("ichorCNA"), sep=".")

    yaml_str <- as.yaml(
        list(
            "${task.process}" = versions
        )
    )
    writeLines(yaml_str, file("versions.yml"))
    """

    stub:
    def args = task.ext.args       ?: ''
    def prefix = task.ext.prefix   ?: "${meta.id}"
    def norm   = normal_wig        ? "normal_wig='${normal_wig}',"          : 'normal_wig=NULL,'
    def pon    = normal_background ? "normal_panel='${normal_background}'," : 'normal_panel=NULL,'
    def map    = map_wig           ? "mapWig='${map_wig}',"                 : 'mapWig=NULL,'
    def centro = centromere        ? "centromere='${centromere}',"          : ''
    def rep    = rep_time_wig      ? "repTimeWig='${rep_time_wig}',"        : 'repTimeWig=NULL,'
    def exons  = exons             ? "exons.bed='${exons}',"                : ''
    """
    #!/usr/bin/env Rscript
    library("ichorCNA")
    library("yaml")

    ### Make Versions YAML for NF-Core ###
    versions = list()
    versions["r"]        <- paste(R.Version()\$major, R.Version()\$minor, sep=".")
    versions["ichorCNA"] <- paste(packageVersion("ichorCNA"), sep=".")

    yaml_str <- as.yaml(
        list(
            "${task.process}" = versions
        )
    )
    writeLines(yaml_str, file("versions.yml"))
    """

}

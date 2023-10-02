process ICHORCNA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-ichorcna=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.5.0--pl5321r42hdfd78af_0' :
        'biocontainers/r-ichorcna:0.5.0--pl5321r42hdfd78af_0' }"

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
    tuple val(meta), path("**.cna.seg")    , emit: cna_seg
    tuple val(meta), path("**.params.txt") , emit: ichorcna_params
    tuple val(meta), path("**.pdf")        , emit: genome_plot
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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
}

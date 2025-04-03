process ICHORCNA_CREATEPON {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.5.1--r43hdfd78af_0' :
        'biocontainers/r-ichorcna:0.5.1--r43hdfd78af_0' }"

    input:
    path wigs
    path gc_wig
    path map_wig
    path centromere
    path rep_time_wig
    path exons

    output:
    path "${prefix}*.rds", emit: rds
    path "${prefix}*.txt", emit: txt
    path "versions.yml"  , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args     ?: ''
    prefix = task.ext.prefix ?: "PoN"
    def map    = map_wig         ? "mapWig='${map_wig}',"                 : 'mapWig=NULL,'
    def centro = centromere      ? "centromere='${centromere}',"          : ''
    def rep    = rep_time_wig    ? "repTimeWig='${rep_time_wig}',"        : 'repTimeWig=NULL,'
    def exons  = exons           ? "exons.bed='${exons}',"                : ''

    """
    #!/usr/bin/env Rscript
    library("ichorCNA")
    library("yaml")

    write.table(strsplit("${wigs}"," ")[[1]],"filelist.txt", row.names = FALSE, col.names = FALSE)

    createPanelOfNormals(
        gcWig='${gc_wig}',
        ${map}
        ${rep}
        filelist = "filelist.txt",
        outfile = "${prefix}",
        ${exons}
        ${centro}
        $args
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
    prefix = task.ext.prefix ?: "PoN"
    """
    #!/usr/bin/env Rscript
    library("yaml")

    file.create("${prefix}.rds")
    file.create("${prefix}.txt")

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

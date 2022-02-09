
process ASCAT {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::ascat=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ascat:3.0.0':
        'quay.io/biocontainers/ascat:3.0.0--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    library(ASCAT)
    ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt","Germline_BAF.txt")
    ascat.plotRawData(ascat.bc)
    ascat.bc = ascat.aspcf(ascat.bc)
    ascat.plotSegmentedData(ascat.bc)
    ascat.output = ascat.runAscat(ascat.bc)


    #version export. Have to hardcode process name and software name because
    #won't run inside an R-block
    version_file_path="versions.yml"
    f <- file(version_file_path,"w")
    writeLines("ASCAT:", f)
    writeLines(paste0(" ascat: 2.5.2",f)
    close(f)

    """
}

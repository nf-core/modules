
process ASCAT {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::ascat=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ascat:3.0.0':
        'biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce' }"

    input:
    tuple val(meta), path(tumorbam)
    tuple val(meta), path(normalbam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    library(ASCAT)
    d<-ascat.prepareHTS(
      tumourseqfile = $tumorbam,
      normalseqfile = $normalbam,
      tumourname = "Tumour",
      normalname = "Normal",
      allelecounter_exe = "alleleCounter",
      alleles.prefix = "G1000_alleles_hg19_chr",
      loci.prefix = "G1000_loci_hg19_chr",
      gender = "XX",
      genomeVersion = "hg19",
      nthreads = 1)





    #version export. Have to hardcode process name and software name because
    #won't run inside an R-block
    version_file_path="versions.yml"
    f <- file(version_file_path,"w")
    writeLines("ASCAT:", f)
    writeLines(paste0(" ascat: 3.0.0",f)
    close(f)

    """
}

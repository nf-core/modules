
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
    tuple val(meta), path("*.rdata"), emit: rdata
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    library(ASCAT)

    #Must download loci and allele file from dropbox and unzip
    #https://www.dropbox.com/s/l3m0yvyca86lpwb/G1000_loci_hg19.zip
    #https://www.dropbox.com/s/3fzvir3uqe3073d/G1000_alleles_hg19.zip


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
      nthreads = $task.cpus
    )


    ascat.bc = ascat.loadData(
      Tumor_LogR_file = "Tumour_tumourLogR.txt",
      Tumor_BAF_file = "Tumour_normalBAF.txt",
      Germline_LogR_file = "Tumour_normalLogR.txt",
      Germline_BAF_file = "Tumour_normalBAF.txt",
      genomeVersion = "hg19"
    )

    ascat.bc = ascat.aspcf(ascat.bc)

    ascat.output = ascat.runAscat(ascat.bc)

    save(ascat.output, file = "ascat.output.rdata)

    #version export. Have to hardcode process name and software name because
    #won't run inside an R-block
    version_file_path="versions.yml"
    f <- file(version_file_path,"w")
    writeLines("ASCAT:", f)
    writeLines(paste0(" ascat: 3.0.0",f)
    close(f)

    """
}

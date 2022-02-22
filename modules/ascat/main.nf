process ASCAT {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::ascat=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ascat:3.0.0':
        'quay.io/biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0' }"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    // path(gcfile) #TODO add in this later, it's a required argument in sarek
    path(allele_files)
    path(loci_files)
    

    output:
    tuple val(meta), path("*.png"),              emit: png
    tuple val(meta), path("*cnvs.txt"),          emit: cnvs
    tuple val(meta), path("*purityploidy.txt"),  emit: purityploidy
    tuple val(meta), path("*segments.txt"),      emit: segments
    path "versions.yml",                         emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def purity = args.purity ? "$args.purity" : "NULL"
    def ploidy = args.ploidy ? "$args.ploidy" : "NULL"
    def gender = args.gender ? "$args.gender" : "'XX'"
    def genomeVersion = args.genomeVersion ? "$args.genomeVersion" : "'hg19'"

    """
    #!/usr/bin/env Rscript

    #TODO: set (temporary) test arguments, figure out how to define later
    gcfile_args = NULL
    ploidy = $ploidy
    purity = $purity
    gender = $gender
    genomeVersion =$genomeVersion

    print("tests")
    print(ploidy)
    print(purity)

    library(RColorBrewer)
    library(ASCAT)
    options(bitmapType='cairo')

    #prepare from BAM files
    ascat.prepareHTS(
      tumourseqfile = "$tumor_bam",
      normalseqfile = "$normal_bam",
      tumourname = "Tumour",
      normalname = "Normal",
      allelecounter_exe = "alleleCounter",
      alleles.prefix = "$allele_files/G1000_alleles_hg19_chr",
      loci.prefix = "$loci_files/G1000_loci_hg19_chr",
      gender = gender,
      genomeVersion = genomeVersion,
      chrom_names = c("21","22"), #TODO: remove this, it's only for testing
      nthreads = $task.cpus
    )

    #Load the data
    ascat.bc = ascat.loadData(
      Tumor_LogR_file = "Tumour_tumourLogR.txt",
      Tumor_BAF_file = "Tumour_normalBAF.txt",
      Germline_LogR_file = "Tumour_normalLogR.txt",
      Germline_BAF_file = "Tumour_normalBAF.txt",
      genomeVersion = genomeVersion
    )

    #GC wave correction (TODO: re-activate - it's a mandatory argument in run_sarek)
    #ascat.bc = ascat.GCcorrect(ascat.bc, gcfile_args)

    #Plot the raw data
    ascat.plotRawData(ascat.bc)

    #Segment the data
    ascat.bc = ascat.aspcf(ascat.bc)

    #Plot the segmented data
    ascat.plotSegmentedData(ascat.bc)

    #Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
    #If psi and rho are manually set:
    if (!is.null(purity) && !is.null(ploidy)){
      ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=purity, psi_manual=ploidy)
    } else if(!is.null(purity) && is.null(ploidy)){
      ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=purity)
    } else if(!is.null(ploidy) && is.null(purity)){
      ascat.output <- ascat.runAscat(ascat.bc, gamma=1, psi_manual=ploidy)
    } else {
      ascat.output <- ascat.runAscat(ascat.bc, gamma=1)
    }

    #Write out segmented regions (including regions with one copy of each allele)
    write.table(ascat.output[["segments"]], file=paste0("$prefix", ".segments.txt"), sep="\t", quote=F, row.names=F)

    #Write out CNVs in bed format
    cnvs=ascat.output[["segments"]][2:6]
    write.table(cnvs, file=paste0("$prefix",".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    #Write out purity and ploidy info
    summary <- tryCatch({
    		matrix(c(ascat.output[["aberrantcellfraction"]], ascat.output[["ploidy"]]), ncol=2, byrow=TRUE)}, error = function(err) {
    			# error handler picks up where error was generated
    			print(paste("Could not find optimal solution:  ",err))
    			return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
    		}
    )
    colnames(summary) <- c("AberrantCellFraction","Ploidy")
    write.table(summary, file=paste0("$prefix",".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    #version export. Have to hardcode process name and software name because
    #won't run inside an R-block
    version_file_path="versions.yml"
    f <- file(version_file_path,"w")
    writeLines("ASCAT:", f)
    writeLines(" ascat: 3.0.0",f)
    close(f)

    """
}

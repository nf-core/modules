process ASCAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ascat=3.0.0 bioconda::cancerit-allelecount-4.3.0": null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0':
        'quay.io/biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0' }"
    
    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    tuple val(meta_allele), path(allele_files)   // needs to be a partially absolute path <path_to_alleles_folder/prefix_without_chr_num_and_ext>
    tuple val(meta_loci), path(loci_files)  

    output:
    tuple val(meta), path("*png"),               emit: png
    tuple val(meta), path("*cnvs.txt"),          emit: cnvs
    tuple val(meta), path("*purityploidy.txt"),  emit: purityploidy
    tuple val(meta), path("*segments.txt"),      emit: segments
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args        ?: ''
    def processName    =  "'" + "${task.process}" + "'"
    def prefix         = task.ext.prefix      ?: "${meta.id}"  
    def allele_prefix  = "\${PWD}/${allele_files}"    
    //print allele_prefix
    def gender         = args.gender          ?  "$args.gender" :        "NULL"
    def genomeVersion  = args.genomeVersion   ?  "$args.genomeVersion" : "NULL"
    def purity         = args.purity          ?  "$args.purity" :        "NULL"
    def ploidy         = args.ploidy          ?  "$args.ploidy" :        "NULL"
    def gc_files       = args.gc_files        ?  "$args.gc_files" :      "NULL"

    def minCounts_arg                    = args.minCounts                     ?  ",minCounts = $args.minCounts" : ""
    def chrom_names_arg                  = args.chrom_names                   ?  ",chrom_names = $args.chrom_names" : ""
    def min_base_qual_arg                = args.min_base_qual                 ?  ",min_base_qual = $args.min_base_qual" : ""
    def min_map_qual_arg                 = args.min_map_qual                  ?  ",min_map_qual = $args.min_map_qual" : ""
    def ref_fasta_arg                    = args.ref_fasta                     ?  ",ref.fasta = '$args.ref_fasta'" : ""
    def skip_allele_counting_tumour_arg  = args.skip_allele_counting_tumour   ?  ",skip_allele_counting_tumour = $args.skip_allele_counting_tumour" : ""
    def skip_allele_counting_normal_arg  = args.skip_allele_counting_normal   ?  ",skip_allele_counting_normal = $args.skip_allele_counting_normal" : ""



    """
    echo ${allele_prefix}
    echo "$allele_prefix"

    
    """


}

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process METAPHLAN3_RUN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::metaphlan=3.0.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/metaphlan:3.0.9--pyhb7b1952_0"
    } else {
        container "quay.io/biocontainers/metaphlan:3.0.9--pyhb7b1952_0"
    }

    input:
    // As metaphlan 3.0 can handle a variety of input types (fastq, fasta, sam and intermediate bt2 map (bowtie2out)) I have created an input_type function to handle this
    tuple val(meta), path(input)
    path metaphlan_db

    output:
    tuple val(meta), path("*_profile.txt")   ,                emit: profile
    tuple val(meta), path("*.biom")          ,                emit: biom
    tuple val(meta), path('*.bowtie2out.txt'), optional:true, emit: bt2out
    path "*.version.txt"                     ,                emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    //contains.() method to handle sam/bam and uncompressed files.. Any better suggestions would be greatly appreciated!!
<<<<<<< HEAD
    def input_type  = ("$input".endsWith(".fastq.gz")) ? "--input_type fastq" :  ("$input".contains(".fasta")) ? "--input_type fasta" : ("$input".endsWith(".bowtie2out.txt")) ? "--input_type bowtie2out" : "--input_type sam" // .contains() to accomodate fasta non-compressed input also?
=======
    def input_type  = ("$input".contains(".fastq")) ? "--input_type fastq" :  ("$input".contains(".fasta")) ? "--input_type fasta" : ("$input".endsWith(".bowtie2out.txt")) ? "--input_type bowtie2out" : "--input_type sam" // .contains() to accomodate non-compressed input?
>>>>>>> ff7e45c8ab3af633135eec76add72d70b3e79134
    def input_data  = ("$input_type".contains("fastq")) && !meta.single_end ? "${input[0]},${input[1]}" : "$input"
    def bowtie2_out = "$input_type" == "--input_type bowtie2out" || "$input_type" == "--input_type sam" ? '' : "--bowtie2out ${prefix}.bowtie2out.txt" // no intermediate alignment files produced for sam & bowtie2out input
    """
    metaphlan \\
        --nproc $task.cpus \\
        $input_type \\
        $input_data \\
        $options.args \\
        $bowtie2_out \\
        --bowtie2db ${metaphlan_db} \\
        --biom ${prefix}.biom \\
        --output_file ${prefix}_profile.txt

    echo \$(metaphlan --version 2>&1) | awk '{print \$3}' > ${software}.version.txt
    """
}

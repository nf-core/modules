/*
========================================================================================
   VerifyBamID2 module
========================================================================================
   Website: http://griffan.github.io/VerifyBamID/
========================================================================================
*/

process VERIFYBAMID_VERIFYBAMID2 {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::verifybamid2=2.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid2:2.0.1--hbb20b25_6' :
        'quay.io/repository/biocontainers/verifybamid2' }"

    input:
    // // Owing to the nature of verifybamid we here provide solutions to working with bam files and optional
    // // alternative variant files, for use with the 'diff' suite of tools.
    // // Other optional input files can be utilised in a similar way to below but we do not exhaustively itterate through all
    // // possible options. Instead we leave that to the user.
    tuple val(meta), path(bam), path(bai)
    // tuple path(svd_prefix.ud), path(svd_prefix.mu), path(svd_bed)
    // path references 

    output:
    // tuple val(meta), path("*.log")                            , emit: log
    // tuple val(meta), path("*.selfSM")          , optional:true, emit: self_SM
    // tuple val(meta), path(“*.Ancestry”)        , optional:true, emit: ancestry
    path "versions.yml"                                       , emit: versions

    // when:
    // task.ext.when == null || task.ext.when

    // script:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // def args_list = args.tokenize()
    
    // def bam_file = ("$bam".endsWith(".bam")) ? "--BamFile ${bam}" : 
    //     ("$bam".endsWith(".cram")) ? "--BamFile ${bam}" : ''
    // args_list.removeIf { it.contains('--BamFile') }

    // def svd_args = "--UDPath ${svd_prefix.svd_ud} --MeanPath ${svd_prefix.svd_mu} --BedPath ${svd_prefix.svd_bed}"

    // if (args.contains('--UDPath') & args.contains('--MeanPath') & args.contains('--BedPath') )
    //     args_list.removeIf { it.contains('--SVDPrefix') }
    //     svd_args = ""

    // def referece_arg = ($args.contains('--Reference')) ? "--Reference ${references}" : 
    //     ("referecens".endsWith(".fasta")) ? "--References ${references}" : ''
    // args_list.removeIf { it.contains('--Reference') }

    // // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    // //               If the software is unable to output a version number on the command-line then it can be manually specified
    // //               e.g. https://github.com/nf-core/modules/blob/master/modules/homer/annotatepeaks/main.nf
    // //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    // //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // // TODO nf-core: Please replace the example samtools command below with your module's command
    // // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    verifybamid2 \\
        // $svd_arg \\
        // $bam_file \\
        // $reference_arg  \\
        // ${args_list.join(' ')} \\
        // > vb_contamination.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid: \$(echo \$(verifybamid2 --help 2>&1) | sed '3,5p;d' | sed 's/ Version://; s/ Copyright/Copyright/; s/ This/This/')
    END_VERSIONS
    """
}

/*
========================================================================================
   VerifyBamID2 module
========================================================================================
   Website: http://griffan.github.io/VerifyBamID/
========================================================================================
*/

process VERIFYBAMID_VERIFYBAMID2 {
    tag '${meta.id}'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::verifybamid2=2.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid2:2.0.1--hbb20b25_6' :
        'quay.io/repository/biocontainers/verifybamid2' }"

    input:
    // Owing to the nature of verifybamid we here provide solutions to working with bam files and optional 
    // alternative variant files. Other optional input files & flags can be utilised in a similar way to below but 
    // we do not exhaustively list all possible options. Instead we leave that to the user.
    tuple val(meta), path(bam), path(bai)
    tuple path(svd_ud), path(svd_mu), path(svd_bed)
    path references 

    output:
    tuple val(meta), path("*.log")             , optional:true, emit: log
    tuple val(meta), path("*.selfSM")          , optional:true, emit: self_SM
    tuple val(meta), path("*.Ancestry")        , optional:true, emit: ancestry
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    
    def bam_file = ("$bam".endsWith('.bam')) ? "--BamFile ${bam}" : 
        ("$bam".endsWith(".cram")) ? "--BamFile ${bam}" : ""
    args_list.removeIf { it.contains('--BamFile') }

    def svd_args = ( svd_ud.baseName.equals(svd_mu.baseName) && svd_ud.baseName.equals(svd_bed.baseName) ) ? 
        "--SVDPrefix ${svd_ud.baseName}" : "${svd_ud.baseName}, ${svd_mu.baseName}, ${svd_bed.baseName}"

    if (args.contains('--UDPath') && args.contains('--MeanPath') && args.contains('--BedPath')) {
        svd_args = svd_args + "--UDPath ${svd_ud} --MeanPath ${svd_mu} --BedPath ${svd_bed}"
    }
    args_list.removeIf { it.contains('--UDPath') }
    args_list.removeIf { it.contains('--MeanPath') }
    args_list.removeIf { it.contains('--BedPath') }

    def reference_args = ("$references".endsWith('.fasta')) ? 
        "--Reference ${references}" : ""

    // Enable generating customized reference stack
    if (args.contains('--RefVCF')) {
        svd_args = ""
        args_list.removeIf { it.contains('--UDPath') }
        args_list.removeIf { it.contains('--MeanPath') }
        args_list.removeIf { it.contains('--BedPath') }
    }

    """
    verifybamid2 \\
        --NumThread $task.cpus \\
        ${svd_args} \\
        ${bam_file} \\
        ${reference_args}  \\
        ${args_list.join(' ')} \\
        > vb_contamination.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid: \$(echo \$(verifybamid2 --help 2>&1 | sed -e '3p;d' | sed -e 's/ Version://'))
    END_VERSIONS
    """
}
process WIPERTOOLS_FASTQSCATTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.5--pyhdfd78af_0':
        'quay.io/biocontainers/wipertools:1.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    val(num_splits)

    output:
    tuple val(meta), path("${out_folder}/*") , emit: fastq_chunks
    tuple val("${task.process}"), val('wipertools'), eval("wipertools fastqscatter --version"), topic: versions, emit: versions_wipertools


    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def args_list   = args.tokenize()
    out_folder      = (args_list.contains('--out_folder') ? args_list[args_list.indexOf('--out_folder')+1] :
                        (args_list.contains('-o') ? args_list[args_list.indexOf('-o')+1] : 'chunks'))
    if(!args.contains('-o') && !args.contains('--out_folder')) {
        args += " -o ${out_folder}"
    }
    """
    wipertools \\
        fastqscatter \\
        -f ${fastq} \\
        -n ${num_splits} \\
        -p ${prefix} \\
        ${args}
    """

    stub:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    out_folder    = (args_list.contains('--out_folder') ? args_list[args_list.indexOf('--out_folder')+1] :
                        (args_list.contains('-o') ? args_list[args_list.indexOf('-o')+1] : 'chunks'))
    if(!args.contains('-o') && !args.contains('--out_folder')) {
        args += " -o ${out_folder}"
    }
    """
    mkdir ${out_folder}
    for i in {1..${num_splits}}
    do
        echo "" | gzip > ${out_folder}/${prefix}_\$i-of-${num_splits}_suffix.fastq.gz
    done
    """
}

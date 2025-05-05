process WIPERTOOLS_FASTQSCATTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.5--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    val(num_splits)

    output:
    tuple val(meta), path("${out_folder}/*") , emit: fastq_chunks
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqscatter: \$(wipertools fastqscatter --version)
    END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqscatter: \$(wipertools fastqscatter --version)
    END_VERSIONS
    """
}

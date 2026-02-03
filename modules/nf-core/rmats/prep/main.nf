process RMATS_PREP {
    tag "${meta.id}"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rmats:4.3.0--py311hf2f0b74_5'
        : 'biocontainers/rmats:4.3.0--py311hf2f0b74_5'}"

    input:
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct

    tuple val(meta), path(genome_bam)
    path reference_gtf
    val rmats_read_len

    output:
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    tuple val(meta), path("${prefix}_prep*.rmats"), emit: prep_rmats_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)

    //   --readLength READLENGTH
    //                    The length of each read. Required parameter, with the
    //                    value set according to the RNA-seq read length
    // TODO - question. Does this definition mean I should change it by read length? If so, look at a samtools command to figure it out. Samtools stats!
    // TODO - should I modify the prefix to include rmats_prep only in a subworkflow via modules.config? It seems so, see example at https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L576
    """
    echo ${genome_bam} > ${prefix}.prep.b1.txt

    rmats \\
        --task prep \\
        ${args} \\
        --nthread ${task.cpus} \\
        -- b1 ${prefix}.prep.b1.txt \\
        --gtf ${reference_gtf} \\
        --readLength ${rmats_read_len} \\
        --tmp ${prefix}_rmats_tmp \\
        --od ${prefix}_rmats_prep

    for file in `ls ${prefix_rmats}_tmp/*.rmats`
    do
    cp ${prefix_rmats}_tmp/${file} ${prefix_}prep_${file}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(rmats.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo ${args}

    touch ${prefix}.rmats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats:  \$(echo \$(rmats.py --version) | sed -e "s/v//g")
    END_VERSIONS
    """
}

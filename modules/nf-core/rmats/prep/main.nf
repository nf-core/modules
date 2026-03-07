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
    // TODO - post seems to need only the BAM *names*, not the actual files. Could we just get the first line of each file to get the names?
    // for file in `ls multi_bam_rmats_prep_tmp/*.rmats`; do head -1 $file; done | tr '\n' ','
    // TODO - for stats, it should be possible to parse the formula using patsy, but if we include PAIRADISE we might have R - just do this in R, first pass
    path reference_gtf
    val rmats_read_len

    output:
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    tuple val(meta), path("*.rmats"), emit: prep_rmats_file
    tuple val(meta), path("*outcomes_by_bam.txt"), emit: prep_read_outcomes_file
    tuple val("${task.process}"), val('rmats'), eval('rmats.py --version | sed -e "s/v//g"'), emit: versions_rmats, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)

    //   --readLength READLENGTH
    //                    The length of each read. Required parameter, with the
    //                    value set according to the RNA-seq read length
    // TODO - question. Does this definition mean I should change it by read length? If so, look at a samtools command to figure it out. Samtools stats!
    // TODO - should I modify the prefix to include rmats_prep only in a subworkflow via modules.config? It seems so, see example at https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L576
    """
    echo ${genome_bam} > ${prefix}.prep.b1.txt

    rmats.py \\
        --task prep \\
        ${args} \\
        --nthread ${task.cpus} \\
        --b1 ${prefix}.prep.b1.txt \\
        --gtf ${reference_gtf} \\
        --readLength ${rmats_read_len} \\
        --tmp ${prefix}_rmats_tmp \\
        --od ${prefix}_rmats_prep

    for file in `ls ${prefix}_rmats_tmp/*`
    do
    cp \${file} ${prefix}_prep_\$(basename \${file})
    done
    """

    // NOTES for post - post requires the rmats files to be in the tmp directory, otherwise it fails

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.rmats
    touch ${prefix}_outcomes_by_bam.txt
    """
}

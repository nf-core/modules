process RMATS_PREP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rmats:4.3.0--py311hf2f0b74_5'
        : 'biocontainers/rmats:4.3.0--py311hf2f0b74_5'}"

    input:
    tuple val(meta), path(genome_bam)
    // NOTES - post seems to need only the BAM *names*, not the actual files. Could we just get the first line of each file to get the names?
    // for file in `ls multi_bam_rmats_prep_tmp/*.rmats`; do head -1 $file; done | tr '\n' ','
    // possible suggestions from @SPPearce - pass ${prefix}.prep.b1.txt as outut
    // NOTES - for stats, it should be possible to parse the formula using patsy, but if we include PAIRADISE we might have R - just do this in R, first pass
    tuple val(meta2), path(reference_gtf)
    val rmats_read_len

    output:
    tuple val(meta), path("*.rmats"), emit: prep_rmats_file
    tuple val(meta), path("*read_outcomes_by_bam.txt"), emit: prep_read_outcomes_file
    tuple val("${task.process}"), val('rmats'), eval('rmats.py --version | sed -e "s/v//g"'), emit: versions_rmats, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // NOTES   --readLength READLENGTH
    //                    The length of each read. Required parameter, with the
    //                    value set according to the RNA-seq read length
    //          I should change it by read length (in workflow)! Look at Samtools stats!
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

    cp ${prefix}_rmats_tmp/*.txt ${prefix}_read_outcomes_by_bam.txt
    cp ${prefix}_rmats_tmp/*.rmats ${prefix}.rmats
    """

    // NOTES for post - post requires the rmats files to be in the tmp directory, otherwise it fails

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.rmats
    touch ${prefix}_read_outcomes_by_bam.txt
    """
}

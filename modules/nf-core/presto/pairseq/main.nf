process PRESTO_PAIRSEQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.9--pyhdfd78af_0':
        'quay.io/biocontainers/presto:0.7.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1_reads), path(R2_reads)
    val(barcode_position)

    output:
    tuple val(meta), path("*_pair-pass.fastq.gz"), path("*_pair-pass.fastq.gz") , emit: reads
    path "*_command_log.txt", emit: logs
    tuple val("${task.process}"), val('presto'), eval('PairSeq.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def copyfield = (barcode_position == 'R1')? '--1f BARCODE' : (barcode_position == 'R2')? '--2f BARCODE' : (barcode_position == 'R1R2')? '--1f BARCODE --2f BARCODE' : (barcode_position == 'clustersets')? '--1f CLUSTER --2f CLUSTER' : ''
    def args = task.ext.args?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    PairSeq.py -1 ${R1_reads} \\
               -2 ${R2_reads} \\
               --outname ${prefix} \\
               $copyfield \\
               $args > ${prefix}_command_log.txt

    """

    stub:
    def copyfield = (barcode_position == 'R1')? '--1f BARCODE' : (barcode_position == 'R2')? '--2f BARCODE' : (barcode_position == 'R1R2')? '--1f BARCODE --2f BARCODE' : (barcode_position == 'clustersets')? '--1f CLUSTER --2f CLUSTER' : ''
    def args = task.ext.args?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-1_pair-pass.fastq.gz \\
          ${prefix}-2_pair-pass.fastq.gz \\
          ${prefix}_command_log.txt 
    """
}

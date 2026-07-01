// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

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
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/624977dfaf562211e68a8a868ca80acc8461f1ac/modules/nf-core/cutadapt/main.nf#L34-L46
    //               Complex example: https://github.com/nf-core/modules/blob/88d43dad73a675e66bff49ebb57fe657a5909018/modules/nf-core/bedtools/split/main.nf#L32-L43
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    touch ${prefix}-1_pair-pass.fastq.gz \\
          ${prefix}-2_pair-pass.fastq.gz \\
          ${prefix}_command_log.txt 
    """
}

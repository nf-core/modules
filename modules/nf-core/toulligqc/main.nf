
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

process TOULLIGQC {
    label 'process_low'

    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/toulligqc:2.5.2--pyhdfd78af_0':
        'biocontainers/toulligqc:2.5.2--pyhdfd78af_0' }"

    input:
    path seq_summary
    path seq_telemetry
    path fast5
    path fastq
    path bam
    
    val report_name
    val output_dir
    path html_report
    path data_report
    path images_dir
    path summary_1dsqr

    val barcodes
    val barcodin
    
    val threads
    val batch_size
    val qscore_threshold
    
    
    output:
    path "*/*.data", emit: report_data
    path "*/*.html", emit: report_html, optional: true
    path "*/images/*.html", emit: plots_html
    path "*/images/plotly.min.js", emit: plotly_js

    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"

    def telemetry_arg = seq_telemetry.name != 'no_telemetry_source' ? "--telemetry-source $seq_telemetry" : ""
    def fast5_arg = fast5.name != 'no_fast5' ? "--fast5-source $fast5" : ""
    def bam_arg = bam.name != 'no_bam' ? "--bam $bam" : ""
    def fastq_arg = fastq.name != 'no_fastq' ? "--fastq $fastq" : ""
    def report_name_arg = report_name != 'no_report_name' ? "-n $report_name" : ""
    def output_arg = output_dir != 'no_output_dir' ? "--output-directory $output_dir" : ""
    def html_arg = html_report.name != 'no_html_report' ? "-o $html_report" : ""
    def data_arg = data_report.name != 'no_data_report' ? "--data-report-path $data_report" : ""
    def images_arg = images_dir.name != 'no_images_dir' ? "--images-directory $images_dir" : ""
    def barcodes_arg = barcodes != 'no_barcodes' ? "--barcodes $barcodes" : ""
    def barcoding_arg = barcodin != 'no_barcoding' ? "--barcoding" : ""
    def threads_arg = threads != 'no_threads' ? "--thread $threads" : ""
    def batch_arg = batch_size != 'no_batch_size' ? "--batch-size $batch_size" : ""
    def qscore_arg = qscore_threshold != 'no_qscore_threshold' ? "--qscore-threshold $qscore_threshold" : ""

    """
    toulligqc -a ${seq_summary}\
        ${telemetry_arg}\
    	${fast5_arg} ${fastq_arg} ${bam_arg}\
    	${report_name_arg}\
    	${output_arg} ${html_arg}\
    	${data_arg} ${images_arg}\
    	${barcoding_arg} ${barcodes_arg}\
    	${threads_arg} ${batch_arg} ${qscore_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        toulligqc: \$(toulligqc --version |& sed '1!d ; s/toulligqc //')
    END_VERSIONS
    """
}

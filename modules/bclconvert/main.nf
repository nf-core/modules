process BCLCONVERT {
    tag 'demultiplexing'
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using bcl-convert. Please use docker or singularity containers."
    }
    container "nfcore/bclconvert:3.9.3"

    input:
    path samplesheet
    path run_dir

    output:
    path "*.fastq.gz"               ,emit: fastq
    path "Reports/*.{csv,xml,bin}"  ,emit: reports
    path "Logs/*.{log,txt}"         ,emit: logs
    path "versions.yml"             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcl-convert \
        $args \\
        --output-directory \$PWD \
        --bcl-input-directory ${run_dir} \\
        --sample-sheet ${samplesheet} \\
        --bcl-num-parallel-tiles $task.cpus \\
        --bcl-num-conversion-threads $task.cpus \\
        --bcl-num-compression-threads $task.cpus \\
        --bcl-num-decompression-threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(echo \$(bcl-convert -V 2>&1) | sed ''s/^.*Version //;s/Copyright.*//'' ))
    END_VERSIONS
    """

    stub:
    """
    touch sample1_S1_L001_R1_001.fastq.gz
    touch sample1_S1_L001_R2_001.fastq.gz
    touch sample1_S1_L002_R1_001.fastq.gz
    touch sample1_S1_L002_R2_001.fastq.gz
    touch sample2_S2_L001_R1_001.fastq.gz
    touch sample2_S2_L001_R2_001.fastq.gz
    touch sample2_S2_L002_R1_001.fastq.gz
    touch sample2_S2_L002_R2_001.fastq.gz

    mkdir Reports
    touch Reports/Adapter_Metrics.csv
    touch Reports/Demultiplex_Stats.csv
    touch Reports/fastq_list.csv
    touch Reports/Index_Hopping_Counts.csv
    touch Reports/IndexMetricsOut.bin
    touch Reports/Quality_Metrics.csv
    touch Reports/RunInfo.xml
    touch Reports/SampleSheet.csv
    touch Reports/Top_Unknown_Barcodes.csv

    mkdir Logs
    touch Logs/Errors.log
    touch Logs/FastqComplete.
    touch Logs/Info.log
    touch Logs/Warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(echo \$(bcl-convert -V 2>&1) | sed ''s/^.*Version //;s/Copyright.*//'' ))
    END_VERSIONS
    """
}

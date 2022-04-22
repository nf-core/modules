process BCLCONVERT {
    tag '$samplesheet'
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
    path "InterOp/*.bin"            ,emit: interop
    path "versions.yml"             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcl-convert \
        $args \\
        --output-directory . \\
        --bcl-input-directory ${run_dir} \\
        --sample-sheet ${samplesheet} \\
        --bcl-num-parallel-tiles ${task.cpus}

    mkdir InterOp
    cp ${run_dir}/InterOp/*.bin InterOp/
    mv Reports/*.bin InterOp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(bcl-convert -V 2>&1 | head -n 1 | sed 's/^.*Version //')
    END_VERSIONS
    """

    stub:
    """
    echo "sample1_S1_L001_R1_001" > sample1_S1_L001_R1_001.fastq.gz
    echo "sample1_S1_L001_R2_001" > sample1_S1_L001_R2_001.fastq.gz
    echo "sample1_S1_L002_R1_001" > sample1_S1_L002_R1_001.fastq.gz
    echo "sample1_S1_L002_R2_001" > sample1_S1_L002_R2_001.fastq.gz
    echo "sample2_S2_L001_R1_001" > sample2_S2_L001_R1_001.fastq.gz
    echo "sample2_S2_L001_R2_001" > sample2_S2_L001_R2_001.fastq.gz
    echo "sample2_S2_L002_R1_001" > sample2_S2_L002_R1_001.fastq.gz
    echo "sample2_S2_L002_R2_001" > sample2_S2_L002_R2_001.fastq.gz

    mkdir Reports
    echo "Adapter_Metrics" >  Reports/Adapter_Metrics.csv
    echo "Demultiplex_Stats" >  Reports/Demultiplex_Stats.csv
    echo "fastq_list" >  Reports/fastq_list.csv
    echo "Index_Hopping_Counts" >  Reports/Index_Hopping_Counts.csv
    echo "IndexMetricsOut" >  Reports/IndexMetricsOut.bin
    echo "Quality_Metrics" >  Reports/Quality_Metrics.csv
    echo "RunInfo" >  Reports/RunInfo.xml
    echo "SampleSheet" >  Reports/SampleSheet.csv
    echo "Top_Unknown_Barcodes" >  Reports/Top_Unknown_Barcodes.csv

    mkdir Logs
    echo "Errors" > Logs/Errors.log
    echo "FastqComplete" > Logs/FastqComplete.txt
    echo "Info" > Logs/Info.log
    echo "Warnings" > Logs/Warnings.log

    mkdir InterOp/
    echo "InterOp" > InterOp/InterOp.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(bcl-convert -V 2>&1 | head -n 1 | sed 's/^.*Version //')
    END_VERSIONS
    """
}

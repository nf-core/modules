process BCLCONVERT {
    tag "${ meta.lane ? meta.id + "." + meta.lane : meta.id }"
    label 'process_high'

    container "nf-core/bclconvert:4.4.6"

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("output/**_S[1-9]*_R?_00?.fastq.gz"), emit: fastq
    tuple val(meta), path("output/**_S[1-9]*_I?_00?.fastq.gz"), emit: fastq_idx, optional: true
    tuple val(meta), path("output/**Undetermined_S0*_R?_00?.fastq.gz"), emit: undetermined, optional: true
    tuple val(meta), path("output/**Undetermined_S0*_I?_00?.fastq.gz"), emit: undetermined_idx, optional: true
    tuple val(meta), path("output/Reports/*.{csv,xml,bin}"), emit: reports
    tuple val(meta), path("output/Logs/*.{log,txt}"), emit: logs
    tuple val(meta), path("output/InterOp/*.bin"), emit: interop, optional: true
    tuple val("${task.process}"), val('bclconvert'), eval("bcl-convert -V 2>&1 | head -n 1 | sed 's/^.*Version //'"), topic: versions, emit: versions_bclconvert

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("BCLCONVERT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def input_tar = run_dir.toString().endsWith(".tar.gz") ? true : false
    def input_dir = input_tar ? run_dir.toString() - '.tar.gz' : run_dir
    """
    if [ ! -d ${input_dir} ]; then
        mkdir -p ${input_dir}
    fi

    if ${input_tar}; then
        ## Ensures --strip-components only applied when top level of tar contents is a directory
        ## If just files or multiple directories, place all in ${input_dir}

        if [[ \$(tar -taf ${run_dir} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
            tar \\
                -C ${input_dir} --strip-components 1 \\
                -xavf \\
                ${args2} \\
                ${run_dir} \\
                ${args3}
        else
            tar \\
                -C ${input_dir} \\
                -xavf \\
                ${args2} \\
                ${run_dir} \\
                ${args3}
        fi
    fi

    bcl-convert \\
        ${args} \\
        --output-directory output \\
        --bcl-input-directory ${input_dir} \\
        --sample-sheet ${samplesheet}

    # copy the InterOp folder contents to ensure it gets picked up when using fusion
    mkdir -p output/InterOp/
    cp -n **/InterOp/*.bin output/InterOp/
    """

    stub:
    """
    mkdir -p output
    echo "fake fastq file" | gzip > output/Sample1_S1_L001_R1_001.fastq.gz
    echo "fake fastq file" | gzip > output/Undetermined_S0_L001_R1_001.fastq.gz

    mkdir -p output/Reports
    echo "fake report file" > output/Reports/Adapter_Cycle_Metrics.csv
    echo "fake report file" > output/Reports/Adapter_Metrics.csv
    echo "fake report file" > output/Reports/Demultiplex_Stats.csv
    echo "fake report file" > output/Reports/Demultiplex_Tile_Stats.csv
    echo "fake report file" > output/Reports/fastq_list.csv
    echo "fake report file" > output/Reports/Index_Hopping_Counts.csv
    echo "fake report file" > output/Reports/IndexMetricsOut.bin
    echo "fake report file" > output/Reports/Quality_Metrics.csv
    echo "fake report file" > output/Reports/Quality_Tile_Metrics.csv
    echo "fake report file" > output/Reports/RunInfo.xml
    echo "fake report file" > output/Reports/SampleSheet.csv
    echo "fake report file" > output/Reports/Top_Unknown_Barcodes.csv

    mkdir -p output/Logs
    echo "fake log file" > output/Logs/Errors.log
    echo "fake log file" > output/Logs/FastqComplete.log
    echo "fake log file" > output/Logs/Info.log
    echo "fake log file" > output/Logs/Warnings.log

    mkdir -p output/InterOp
    echo "fake InterOp file" > output/InterOp/ControlMetricsOut.bin
    echo "fake InterOp file" > output/InterOp/CorrectedIntMetricsOut.bin
    echo "fake InterOp file" > output/InterOp/ErrorMetricsOut.bin
    echo "fake InterOp file" > output/InterOp/ExtractionMetricsOut.bin
    echo "fake InterOp file" > output/InterOp/IndexMetricsOut.bin
    echo "fake InterOp file" > output/InterOp/QMetricsOut.bin
    echo "fake InterOp file" > output/InterOp/TileMetricsOut.bin
    """
}

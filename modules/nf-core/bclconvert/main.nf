process BCLCONVERT {
    tag {"$meta.lane" ? "$meta.id"+"."+"$meta.lane" : "$meta.id" }
    label 'process_high'

    container "nf-core/bclconvert:4.3.13"

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("output/**_S[1-9]*_R?_00?.fastq.gz")        , emit: fastq
    tuple val(meta), path("output/**_S[1-9]*_I?_00?.fastq.gz")        , emit: fastq_idx       , optional:true
    tuple val(meta), path("output/**Undetermined_S0*_R?_00?.fastq.gz"), emit: undetermined    , optional:true
    tuple val(meta), path("output/**Undetermined_S0*_I?_00?.fastq.gz"), emit: undetermined_idx, optional:true
    tuple val(meta), path("output/Reports")                           , emit: reports
    tuple val(meta), path("output/Logs")                              , emit: logs
    tuple val(meta), path("output/InterOp/*.bin")                     , emit: interop         , optional:true
    path("versions.yml")                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BCLCONVERT module does not support Conda. Please use Docker / Singularity / Podman instead."
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
        ## If just files or multiple directories, place all in $input_dir

        if [[ \$(tar -taf ${run_dir} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
            tar \\
                -C $input_dir --strip-components 1 \\
                -xavf \\
                $args2 \\
                $run_dir \\
                $args3
        else
            tar \\
                -C $input_dir \\
                -xavf \\
                $args2 \\
                $run_dir \\
                $args3
        fi
    fi

    bcl-convert \\
        $args \\
        --output-directory output \\
        --bcl-input-directory ${input_dir} \\
        --sample-sheet ${samplesheet}

    # copy the InterOp folder contents to ensure it gets picked up when using fusion
    mkdir -p output/InterOp/
    cp -n **/InterOp/*.bin output/InterOp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(bcl-convert -V 2>&1 | head -n 1 | sed 's/^.*Version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/Reports
    mkdir -p output/Logs
    echo "" | gzip > output/Sample1_S1_L001_R1_001.fastq.gz
    echo "" | gzip > output/Undetermined_S0_L001_R1_001.fastq.gz
    touch output/Reports/Adapter_Cycle_Metrics.csv
    touch output/Reports/Adapter_Metrics.csv
    touch output/Reports/Demultiplex_Stats.csv
    touch output/Reports/Demultiplex_Tile_Stats.csv
    touch output/Reports/fastq_list.csv
    touch output/Reports/Index_Hopping_Counts.csv
    touch output/Reports/IndexMetricsOut.bin
    touch output/Reports/Quality_Metrics.csv
    touch output/Reports/Quality_Tile_Metrics.csv
    touch output/Reports/RunInfo.xml
    touch output/Reports/SampleSheet.csv
    touch output/Reports/Top_Unknown_Barcodes.csv
    touch output/Logs/Errors.log
    touch output/Logs/FastqComplete.log
    touch output/Logs/Info.log
    touch output/Logs/Warnings.log
    mkdir -p output/InterOp
    touch output/InterOp/ControlMetricsOut.bin
    touch output/InterOp/CorrectedIntMetricsOut.bin
    touch output/InterOp/ErrorMetricsOut.bin
    touch output/InterOp/ExtractionMetricsOut.bin
    touch output/InterOp/IndexMetricsOut.bin
    touch output/InterOp/QMetricsOut.bin
    touch output/InterOp/TileMetricsOut.bin
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(bcl-convert -V 2>&1 | head -n 1 | sed 's/^.*Version //')
    END_VERSIONS
    """

}

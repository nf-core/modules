process MGIKIT_DEMULTIPLEX {
    tag {"$run_id"}
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgikit:1.0.0--h3ab6199_0' :
        'biocontainers/mgikit:1.0.0--h3ab6199_0' }"

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("${prefix}/*.fastq.gz")                                                    , emit: fastq
    tuple val(meta), path("${prefix}_undetermined/*.fastq.gz")                                       , optional:true, emit: undetermined
    tuple val(meta), path("${prefix}_ambiguous/*.fastq.gz")                                          , optional:true, emit: ambiguous
    tuple val(meta), path("${prefix}/*mgikit.undetermined_barcode*")                                 , emit: undetermined_reports, optional:true
    tuple val(meta), path("${prefix}/*mgikit.ambiguous_barcode*")                                    , emit: ambiguous_reports, optional:true
    tuple val(meta), path("${prefix}/*mgikit.general")                                               , emit: general_info_reports
    tuple val(meta), path("${prefix}/*mgikit.info")                                                  , emit: index_reports
    tuple val(meta), path("${prefix}/*mgikit.sample_stats")                                          , emit: sample_stat_reports
    tuple val(meta), path("${prefix}/*mgikit.{info,general,ambiguous_barcode,undetermined_barcode}") , emit: qc_reports
    path("versions.yml")                                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    run_id = meta.lane ? "${meta.id}-${meta.lane}" : "${meta.id}"
    prefix = task.ext.prefix ?: "out-${run_id}"

    """
    mgikit demultiplex \\
        -i "${run_dir}" \\
        -s "${samplesheet}" \\
        -o "${prefix}" \\
        ${args}

    if find ${prefix} -name 'Undetermined*.fastq.gz' -print -quit | grep -q .; then
        mkdir -p "${prefix}_undetermined"
        mv ${prefix}/Undetermined*.fastq.gz ${prefix}_undetermined/
    fi

    if find ${prefix} -name 'Ambiguous*.fastq.gz' -print -quit | grep -q .; then
        mkdir -p "${prefix}_ambiguous"
        mv ${prefix}/Ambiguous*.fastq.gz ${prefix}_ambiguous/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgikit: \$(mgikit --version 2>&1 | grep 'MGIKIT - MGI data demultipexing kit' | sed -e 's/MGIKIT - MGI data demultipexing kit. //g')
    END_VERSIONS
    """

    stub:
    run_id = meta.lane ? "${meta.id}-${meta.lane}" : "${meta.id}"
    prefix = task.ext.prefix ?: "out-${run_id}"
    """
    mkdir "${prefix}"
    mkdir -p "${prefix}_undetermined"

    touch "${prefix}/FC1.L01.mgikit.general"
    touch "${prefix}/FC1.L01.mgikit.info"
    touch "${prefix}/FC1.L01.mgikit.undetermined_barcode"
    touch "${prefix}/FC1.L01.mgikit.sample_stats"

    echo "@R001:0001:FC1:1:60:1:3 1:N:0:GACGAATG\\nNNNNNNNN\\n+\\nDDDDDDDD" | gzip > "${prefix}/23-001_S1_L01_R1_001.fastq.gz"
    echo "@R001:0001:FC1:1:60:1:3 2:N:0:GACGAATG\\nNNNNNNNN\\n+\\nDDDDDDDD" | gzip > "${prefix}/23-001_S1_L01_R2_001.fastq.gz"
    echo "@R001:0001:FC1:1:60:1:3 1:N:0:GACGAATG\\nNNNNNNNN\\n+\\nDDDDDDDD" | gzip > "${prefix}/23-002_S2_L01_R1_001.fastq.gz"
    echo "@R001:0001:FC1:1:60:1:3 2:N:0:GACGAATG\\nNNNNNNNN\\n+\\nDDDDDDDD" | gzip > "${prefix}/23-002_S2_L01_R2_001.fastq.gz"

    echo "@R001:0001:FC1:1:60:1:3 1:N:0:GACGAATG\\nNNNNNNNN\\n+\\nDDDDDDDD" | gzip > "${prefix}_undetermined/Undetermined_L01_R1_001.fastq.gz"
    echo "@R001:0001:FC1:1:60:1:3 2:N:0:GACGAATG\\nNNNNNNNN\\n+\\nDDDDDDDD" | gzip > "${prefix}_undetermined/Undetermined_L01_R2_001.fastq.gz"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgikit: \$(mgikit --version 2>&1 | grep 'MGIKIT - MGI data demultipexing kit' | sed -e 's/MGIKIT - MGI data demultipexing kit. //g')
    END_VERSIONS
    """
}

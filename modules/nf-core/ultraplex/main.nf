process ULTRAPLEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ultraplex:1.2.5--py38h4a8c8d9_0' :
        'biocontainers/ultraplex:1.2.5--py38h4a8c8d9_0' }"

    input:
    tuple val(meta), path(reads)
    path(barcode_file)
    val(adapter_seq)

    output:
    tuple val(meta), path("*_matched.fastq.gz")    , emit: fastq
    tuple val(meta), path("*_no_match_*.fastq.gz") , emit: no_match_fastq, optional: true
    path "*.log"                                   , emit: report
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = "1.2.5" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    def adapter_seq_command = ''
    if(adapter_seq) {
        adapter_seq_command = "--adapter ${adapter_seq}"
    }

    read_list = reads.collect{it.toString()}
    if (read_list.size > 1){
        ultraplex_command = """ultraplex \\
        --inputfastq ${read_list[0]} \\
        --input_2 ${read_list[1]} \\
        --barcodes $barcode_file \\
        --threads $task.cpus \\
        --outputprefix ${prefix} \\
        $args \\
        $adapter_seq_command"""
    } else {
        ultraplex_command = """ultraplex \\
        --inputfastq ${read_list[0]} \\
        --barcodes $barcode_file \\
        --threads $task.cpus \\
        --outputprefix ${prefix} \\
        $args \\
        $adapter_seq_command"""
    }

    """
    ${ultraplex_command}

    ## rename the matched files to be able to emit only the matched files
    MATCHES=\$( find . -type f -name '*.fastq.gz' ! -name '*_no_match_*.fastq.gz' )
    for MATCH in \$MATCHES; do
        mv \$MATCH \${MATCH/.fastq.gz/_matched.fastq.gz}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultraplex: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "1.2.5" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ultraplex_${prefix}_Sample1_Fwd_matched.fastq.gz
    echo "" | gzip > ultraplex_${prefix}_Sample1_Rev_matched.fastq.gz
    echo "" | gzip > ultraplex_${prefix}_Sample2_Fwd_matched.fastq.gz
    echo "" | gzip > ultraplex_${prefix}_Sample2_Rev_matched.fastq.gz
    echo "" | gzip > ultraplex_${prefix}_no_match_Fwd.fastq.gz
    echo "" | gzip > ultraplex_${prefix}_no_match_Rev.fastq.gz

    touch ultraplex.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultraplex: $VERSION
    END_VERSIONS
    """
}

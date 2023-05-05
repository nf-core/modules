process ARCASHLA_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::arcas-hla=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.5.0--hdfd78af_0':
        'biocontainers/arcas-hla:0.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fq.gz")                 , emit: extracted_reads_fastq
    path "*.log"                                     , emit: log
    tuple val(meta), path("temp_files/**.sam")       , emit: intermediate_sam,        optional: true
    tuple val(meta), path("temp_files/**.bam")       , emit: intermediate_bam,        optional: true
    tuple val(meta), path("temp_files/**.sorted.bam"), emit: intermediate_sorted_bam, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single_end  = meta.single_end ? "--single" : ""
    def VERSION = "0.5.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    arcasHLA \\
        extract \\
        $args \\
        -t $task.cpus \\
        -o . \\
        --temp temp_files/ \\
        --log ${prefix}.log \\
        $single_end \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: $VERSION
    END_VERSIONS
    """
}

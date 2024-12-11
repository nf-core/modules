process MAGECK_COUNT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mageck=0.5.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9.5--py39h1f90b4d_3':
        'biocontainers/mageck:0.5.9.5--py39h1f90b4d_3' }"

    input:
    tuple val(meta), path(inputfile)
    path(library)

    output:
    tuple val(meta), path("*count.txt"),             emit: count
    tuple val(meta), path("*.count_normalized.txt"), emit: norm
    tuple val(meta), path("*.countsummary.txt"),     emit: summary
    tuple val(meta), path("*.count_normalized.txt"), emit: normalized
    tuple val(meta), path("*.log"),                  emit: logs
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_label = ("$inputfile".endsWith(".fastq.gz") || "$inputfile".endsWith(".fq.gz")) ? "--sample-label ${meta.id}" : ''

    if (meta.single_end && ("$inputfile".endsWith(".fastq.gz") || "$inputfile".endsWith(".fq.gz"))) {
        input = "--fastq ${inputfile}"
    } else {
        input = "--fastq ${inputfile[0]} --fastq-2 ${inputfile[1]}"
    }

    """
    mageck \\
        count \\
        $args \\
        -l $library \\
        -n $prefix \\
        $sample_label \\
        $input


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_file = ("$inputfile".endsWith(".fastq.gz")) ? "--fastq ${inputfile}" : "-k ${inputfile}"
    def sample_label = ("$inputfile".endsWith(".fastq.gz") || "$inputfile".endsWith(".fq.gz")) ? "--sample-label ${meta.id}" : ''
    """
    touch ${prefix}.count.txt
    touch ${prefix}.count_normalized.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}

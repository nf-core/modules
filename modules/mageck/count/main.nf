process MAGECK_COUNT {
    tag "$meta.id"
    label 'process_medium'


    conda (params.enable_conda ? "bioconda::mageck=0.5.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9--py37h6bb024c_0':
        'quay.io/biocontainers/mageck:0.5.9--py37h6bb024c_0' }"

    input:
    tuple val(meta), path(inptfile)
    path(library)
    val(name)

    output:
    tuple val(meta), path("*count*.txt"), emit: count
    path "versions.yml"           , emit: versions
    path("*.count_normalized.txt"), emit: norm

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = [".txt", ".csv", ".tsv"]
    def input_file = ("$inptfile".endsWith(".fastq.gz")) ? "--fastq ${inptfile}" :
        ("$inptfile".endsWith(tuple(ext)) ? "-k ${inptfile}" : ''
    def sample_label = ("$inptfile".endsWith(".fastq.gz")) ? "--sample-label ${meta.id}" : ''
    """

    mageck \\
        count \\
        $args \\
        -l $library \\
        -n $name \\
        $sample_label \\
        $input_file \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}

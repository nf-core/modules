
process LIGHTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lighter:1.1.3--hdcf5f25_0':
        'biocontainers/lighter:1.1.3--hdcf5f25_0' }"

    input:
    tuple val(meta), path(fastq_folder, stageAs: "?/*.fastq.gz") 
    val genome_size //Estimated or user specified genome size
    val kmer_size //K_mer size
    val alpha // define alpha is not necessary, because "When using "-K" instead of "-k", Lighter will go through the reads an extra pass to decide C. "
    val output_dir

    output:
    tuple val(meta), path("*.fq.gz"), emit: outfastq
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def k_par = $alpha ? "-k ${kmer_size} ${genome_size} ${alpha}" : "-K ${kmer_size} ${genome_size}"
    def input_fastq = ${fastq_folder.listFiles().collect { "-r $it" }.join(' ')}

    """
    lighter \\
        ${input_fastq} \\
        $k_par \\
        -t $task.cpus \\
        -od ${output_dir} \\
        $args \\
     

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lighter: \$(lighter -v |& sed '1!d ; s/lighter //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  
    """
    touch ${prefix}.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lighter: \$(lighter -v |& sed '1!d ; s/lighter //')
    END_VERSIONS
    """
}

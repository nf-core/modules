process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim*.fastq.gz")          , emit: trimmed_reads
    tuple val(meta), path("*.unpaired.trim_*.fastq.gz"), optional:true, emit: unpaired_reads
    tuple val(meta), path("*.log")                     , emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "SE" : "PE"
    def output = meta.single_end ?
        "${prefix}.trim.fastq.gz"
        : "${prefix}.paired.trim_1.fastq.gz ${prefix}.unpaired.trim_1.fastq.gz ${prefix}.paired.trim_2.fastq.gz ${prefix}.unpaired.trim_2.fastq.gz"
    // FIXME Handle SE defaults
    def qual_trim = task.ext.args2 ?: "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    """
    trimmomatic \\
        $trimmed \\
        -threads $task.cpus \\
        -trimlog ${prefix}.log \\
        $reads \\
        $output \\
        $qual_trim \\
        $args
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}

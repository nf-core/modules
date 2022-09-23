process GATK4_COLLECTREADCOUNTS {
    tag "$meta.id"
    label 'process_medium'

<<<<<<< HEAD
    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.0--hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.2.6.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    path  intervals
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path  "versions.yml"          , emit: versions
=======
    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("*.hdf5"), optional: true, emit: hdf5
    tuple val(meta), path("*.tsv") , optional: true, emit: tsv
    path "versions.yml"           , emit: versions
>>>>>>> master

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
<<<<<<< HEAD
    def interval_command = intervals ? "--intervals $intervals" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK CollectReadCounts] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
=======

    def reference = fasta ? "--reference $fasta" : ""
    def extension = args.contains("--format HDF5") ? "hdf5" :
                    args.contains("--format TSV")  ? "tsv" :
                    "hdf5"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK COLLECTREADCOUNTS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
>>>>>>> master
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" CollectReadCounts \\
        --input $input \\
<<<<<<< HEAD
        --reference $fasta \\
        $interval_command \\
        -imr OVERLAPPING_ONLY \\
        --format TSV \\
        -O ${meta.id}.tsv \\
=======
        --intervals $intervals \\
        --output ${prefix}.$extension \\
        $reference \\
        --tmp-dir . \\
>>>>>>> master
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

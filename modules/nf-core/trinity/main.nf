process TRINITY {
    tag "$meta.id"
    label 'process_high_memory'

    conda "bioconda::trinity=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
        'biocontainers/trinity:2.13.2--h00214ad_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fa.gz")       , emit: transcript_fasta
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        reads_args = "--single ${reads}"
    } else {
        reads_args = "--left ${reads[0]} --right ${reads[1]}"
    }

    // --seqType argument, fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    seqType_args = reads[0] ==~ /(.*fasta(.gz)?$)|(.*fa(.gz)?$)/ ? "fa" : "fq"

    // Define the memory requirements. Trinity needs this as an option.
    def avail_mem = 7
    if (!task.memory) {
        log.info '[Trinity] Available memory not known - defaulting to 7GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

    """
    # Note that Trinity needs the word 'trinity' in the outdir

    Trinity \\
    --seqType ${seqType_args} \\
    --max_memory ${avail_mem}G \\
    ${reads_args} \\
    --output ${prefix}_trinity \\
    --CPU $task.cpus \\
    $args

    gzip -cf ${prefix}_trinity.Trinity.fasta > ${prefix}.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: \$(echo \$(Trinity --version | head -n 1 2>&1) | sed 's/^Trinity version: Trinity-v//' ))
    END_VERSIONS

    # Need to only take the first line of --version since it will warn about not being up-to-date and this messes up the version.yaml.
    """
}

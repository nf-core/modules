process LIMA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::lima=2.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lima:2.2.0--h9ee0642_0' :
        'quay.io/biocontainers/lima:2.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(ccs)
    path primers

    output:
    tuple val(meta), path("*.counts") , emit: counts
    tuple val(meta), path("*.report") , emit: report
    tuple val(meta), path("*.summary"), emit: summary
    path "versions.yml"               , emit: versions

    tuple val(meta), path("*.bam")              , optional: true, emit: bam
    tuple val(meta), path("*.bam.pbi")          , optional: true, emit: pbi
    tuple val(meta), path("*.{fa, fasta}")      , optional: true, emit: fasta
    tuple val(meta), path("*.{fa.gz, fasta.gz}"), optional: true, emit: fastagz
    tuple val(meta), path("*.fastq")            , optional: true, emit: fastq
    tuple val(meta), path("*.fastq.gz")         , optional: true, emit: fastqgz
    tuple val(meta), path("*.xml")              , optional: true, emit: xml
    tuple val(meta), path("*.json")             , optional: true, emit: json
    tuple val(meta), path("*.clips")            , optional: true, emit: clips
    tuple val(meta), path("*.guess")            , optional: true, emit: guess

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$ccs" == "${prefix}.bam" )      error "Input and output names are the same, set prefix in module configuration"
    if( "$ccs" == "${prefix}.fasta" )    error "Input and output names are the same, set prefix in module configuration"
    if( "$ccs" == "${prefix}.fasta.gz" ) error "Input and output names are the same, set prefix in module configuration"
    if( "$ccs" == "${prefix}.fastq" )    error "Input and output names are the same, set prefix in module configuration"
    if( "$ccs" == "${prefix}.fastq.gz" ) error "Input and output names are the same, set prefix in module configuration"

    """
    OUT_EXT=""

    if [[ $ccs =~ bam\$ ]]; then
        OUT_EXT="bam"
    elif [[ $ccs =~ fasta\$ ]]; then
        OUT_EXT="fasta"
    elif [[ $ccs =~ fasta.gz\$ ]]; then
        OUT_EXT="fasta.gz"
    elif [[ $ccs =~ fastq\$ ]]; then
        OUT_EXT="fastq"
    elif [[ $ccs =~ fastq.gz\$ ]]; then
        OUT_EXT="fastq.gz"
    fi

    lima \\
        $ccs \\
        $primers \\
        $prefix.\$OUT_EXT \\
        -j $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lima: \$( lima --version | sed 's/lima //g' | sed 's/ (.\\+//g' )
    END_VERSIONS
    """
}

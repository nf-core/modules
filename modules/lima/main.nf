process LIMA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::lima=2.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/lima:2.2.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/lima:2.2.0--h9ee0642_0"
    }

    input:
    tuple val(meta), path(ccs)
    path primers

    output:
    tuple val(meta), path("*.clips")  , emit: clips
    tuple val(meta), path("*.counts") , emit: counts
    tuple val(meta), path("*.guess")  , emit: guess
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

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
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

    echo \$OUT_EXT
    lima \\
        $ccs \\
        $primers \\
        $prefix.\$OUT_EXT \\
        -j $task.cpus \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( lima --version | sed 's/lima //g' | sed 's/ (.\\+//g' )
    END_VERSIONS
    """
}

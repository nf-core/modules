process QUALIMAP_BAMQC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::qualimap=2.2.2d" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1"
    } else {
        container "quay.io/biocontainers/qualimap:2.2.2d--1"
    }

    input:
    tuple val(meta), path(bam)
    path gff
    val use_gff

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "versions.yml"              , emit: versions

    script:
    prefix         = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'
    def memory     = task.memory.toGiga() + "G"
    def regions = use_gff ? "--gff $gff" : ''

    def strandedness = 'non-strand-specific'
    if (meta.strandedness == 'forward') {
        strandedness = 'strand-specific-forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'strand-specific-reverse'
    }
    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=$memory \\
        bamqc \\
        $args \\
        -bam $bam \\
        $regions \\
        -p $strandedness \\
        $collect_pairs \\
        -outdir $prefix \\
        -nt $task.cpus

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """
}

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QUALIMAP_BAMQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::qualimap=2.2.2d" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1"
    } else {
        container "quay.io/biocontainers/qualimap:2.2.2d--1"
    }

    input:
    tuple val(meta), path(bam)
    path regions

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "*.version.txt"             , emit: version

    script:
    def software   = getSoftwareName(task.process)
    prefix         = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def memory     = task.memory.toGiga() + "G"
    def threads    = task.cpus
    def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'

    // currently there is no convenient way to have optional input in dsl2
    // there needs to be a placeholder file in the workflow definition with the magic name 'dummy_file.txt' 
    def regions_file = regions.name != 'dummy_file.txt' ? '--gff ${regions}' : ''

    def gcref = ''
    if (meta.genome.toString().startsWith('GRCh')) {
        gcref = '-gd HUMAN'
    } else if (meta.genome.toString().startsWith('GRCm')) {
        gcref = '-gd MOUSE'
    }

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
        $options.args \\
        -bam $bam \\
        $regions_file \\
        -p $strandedness \\
        $collect_pairs \\
        -outdir $prefix \\
        -nt $threads

    echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//' > ${software}.version.txt
    """
}

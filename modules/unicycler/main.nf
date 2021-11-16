// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNICYCLER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::unicycler=0.4.8' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/unicycler:0.4.8--py38h8162308_3"
    } else {
        container "quay.io/biocontainers/unicycler:0.4.8--py38h8162308_3"
    }

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path('*.scaffolds.fa.gz'), emit: scaffolds
    tuple val(meta), path('*.assembly.gfa.gz'), emit: gfa
    tuple val(meta), path('*.log')            , emit: log
    path  "versions.yml"                      , emit: versions

    script:
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def short_reads = shortreads ? ( meta.single_end ? "-s $shortreads" : "-1 ${shortreads[0]} -2 ${shortreads[1]}" ) : ""
    def long_reads  = longreads ? "-l $longreads" : ""
    """
    unicycler \\
        --threads $task.cpus \\
        $options.args \\
        $short_reads \\
        $long_reads \\
        --out ./

    mv assembly.fasta ${prefix}.scaffolds.fa
    gzip -n ${prefix}.scaffolds.fa
    mv assembly.gfa ${prefix}.assembly.gfa
    gzip -n ${prefix}.assembly.gfa
    mv unicycler.log ${prefix}.unicycler.log

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//')
    END_VERSIONS
    """
}

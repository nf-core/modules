// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GUNC_RUN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gunc=1.0.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gunc:1.0.4--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/gunc:1.0.4--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("*maxCSS_level.tsv")                  , emit: maxcss_level_tsv
    tuple val(meta), path("*all_levels.tsv")    , optional: true, emit: all_levels_tsv

    path "versions.yml"                                         , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gunc \\
        run \\
        --input_fasta $fasta \\
        --db_file $db \\
        --threads $task.cpus \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( gunc --version )
    END_VERSIONS
    """
}

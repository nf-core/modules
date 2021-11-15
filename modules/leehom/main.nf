// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION="1.2.15"

process LEEHOM {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::leehom=1.2.15" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/leehom:1.2.15--h29e30f7_1"
    } else {
        container "quay.io/biocontainers/leehom:1.2.15--h29e30f7_1"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.bam")          , optional: true, emit: bam
    tuple val(meta), path("${prefix}.fq.gz")        , optional: true, emit: fq_pass
    tuple val(meta), path("${prefix}.fail.fq.gz")   , optional: true, emit: fq_fail
    tuple val(meta), path("${prefix}_r1.fq.gz")     , optional: true, emit: unmerged_r1_fq_pass
    tuple val(meta), path("${prefix}_r1.fail.fq.gz"), optional: true, emit: unmerged_r1_fq_fail
    tuple val(meta), path("${prefix}_r2.fq.gz")     , optional: true, emit: unmerged_r2_fq_pass
    tuple val(meta), path("${prefix}_r2.fail.fq.gz"), optional: true, emit: unmerged_r2_fq_fail
    tuple val(meta), path("*.log")                                  , emit: log

    path "versions.yml"                                             , emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if ( reads.toString().endsWith('.bam') ) {
            """
            leeHom \\
                $options.args \\
                -t $task.cpus \\
                -o ${prefix}.bam \\
                --log ${prefix}.log \\
                $reads

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                ${getSoftwareName(task.process)}: \$( echo $VERSION )
            END_VERSIONS
            """
    } else if ( meta.single_end ) {
            """
            leeHom \\
                $options.args \\
                -t $task.cpus \\
                -fq1 $reads \\
                -fqo ${prefix} \\
                --log ${prefix}.log

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                ${getSoftwareName(task.process)}: \$( echo $VERSION )
            END_VERSIONS
            """
    } else {
            """
            leeHom \\
                $options.args \\
                -t $task.cpus \\
                -fq1 ${reads[0]} \\
                -fq2 ${reads[1]} \\
                -fqo ${prefix} \\
                --log ${prefix}.log

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                ${getSoftwareName(task.process)}: \$( echo $VERSION )
            END_VERSIONS
            """
    }
}

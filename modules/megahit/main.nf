// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEGAHIT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::megahit=1.2.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/megahit:1.2.9--h2e03b76_1"
    } else {
        container "quay.io/biocontainers/megahit:1.2.9--h2e03b76_1"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("megahit_out/*.contigs.fa")                               , emit: contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.contigs.fa")         , emit: kcontigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.addi.fa")            , emit: addi_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.local.fa")           , emit: localcontigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.final.contigs.fa")   , emit: kfinalcontigs
    path "versions.yml"                                                             , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (meta.single_end) {
        """
        megahit \\
            -r ${reads[0]} \\
            -t $task.cpus \\
            $options.args \\
            --out-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                ${getSoftwareName(task.process)}: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    } else {
        """
        megahit \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -t $task.cpus \\
            $options.args \\
            --out-prefix ${prefix}

        cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                ${getSoftwareName(task.process)}: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
}

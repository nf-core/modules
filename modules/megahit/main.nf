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

    conda (params.enable_conda ? "bioconda::megahit=1.2.9 conda-forge::pigz=2.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("megahit_out/*.contigs.fa.gz")                            , emit: contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    path "versions.yml"                                                             , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        megahit \\
            -r ${reads} \\
            -t $task.cpus \\
            $options.args \\
            --out-prefix $prefix

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $options.args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa

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
            --out-prefix $prefix

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $options.args2 \\
            megahit_out/*.fa \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
}

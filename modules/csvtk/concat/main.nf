// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
params.options = [:]
options        = initOptions(params.options)

process CSVTK_CONCAT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::csvtk=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0"
    }

    input:
    tuple val(meta), path(csv)
    val(in_format)
    val(out_format)

    output:
    tuple val(meta), path("*.${out_format}"), emit: table
    path "versions.yml"                     , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def delimiter = in_format == "tsv" ? "\t" : ","
    def out_delimiter = out_format == "tsv" ? "\t" : ","
    """
    csvtk \\
        concat \\
        $options.args \\
        --num-cpus ${task.cpus} \\
        --delimiter "${delimiter}" \\
        --out-delimiter "${out_delimiter}" \\
        --out-file ${prefix}.${out_format} \\
        ${csv}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}

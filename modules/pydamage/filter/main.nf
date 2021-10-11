// Import generic module functions
<<<<<<< HEAD
include { initOptions; saveFiles; getSoftwareName } from './functions'
=======
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
>>>>>>> c19671dca974354978c9bc1711fca6fe681bdb0b

params.options = [:]
options        = initOptions(params.options)

process PYDAMAGE_FILTER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pydamage=0.62" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pydamage:0.62--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/pydamage:0.62--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("pydamage_results/pydamage_filtered_results.csv"), emit: csv
<<<<<<< HEAD
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
=======
    path "versions.yml"           , emit: versions

    script:
>>>>>>> c19671dca974354978c9bc1711fca6fe681bdb0b
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """

    pydamage \\
        filter \\
        $options.args \\
        $csv

<<<<<<< HEAD
    echo \$(pydamage --version 2>&1)  | sed -e 's/pydamage, version //g' > ${software}.version.txt
=======
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(pydamage --version 2>&1) | sed -e 's/pydamage, version //g')
    END_VERSIONS
>>>>>>> c19671dca974354978c9bc1711fca6fe681bdb0b
    """
}

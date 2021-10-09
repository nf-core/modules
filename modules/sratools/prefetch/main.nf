// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SRATOOLS_PREFETCH {
    tag "$id"

    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0'
    } else {
        container 'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0'
    }

    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), path("$id"), emit: sra
    path "versions.yml"         , emit: versions

    script:
    def software = getSoftwareName(task.process)
    """
    ncbi_settings=\$(vdb-config -o n NCBI_SETTINGS)
    if [[ ! -f "\${ncbi_settings}" ]]; then
        mkdir -p "\${ncbi_settings}"
        echo "${options.vdb_config}" > "\${ncbi_settings}"
    fi

    prefetch \\
        $options.args \\
        --progress \\
        $id

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( prefetch --version 2>&1 | sed 's/^.*version //; s/.*\$//' )
    END_VERSIONS
    """
}

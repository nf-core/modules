// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)
options.vdb_config = params.options.vdb_config ?: "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"

process SRATOOLS_FASTERQDUMP {
    tag "$meta.id"
    label 'process_medium'

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
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path("*.fastq"), emit: reads
    path "versions.yml"             , emit: versions

    script:
    /***************************************************************************
     * N.B.: `fasterq-dump` works fastest when using a memory-mapped directory
     * as temporary workspace. `/tmp` is typically such a directory which is why
     * we explicitly choose it here. You may want to configure `docker.temp` and
     * set it to a memory mapped directory available on the host system.
     **************************************************************************/
    //

    """
    eval "\$(vdb-config -o n NCBI_SETTINGS | sed 's/[" ]//g')"
    if [[ ! -f "\${NCBI_SETTINGS}" ]]; then
        mkdir -p "\$(dirname "\${NCBI_SETTINGS}")"
        printf '${options.vdb_config}' > "\${NCBI_SETTINGS}"
    fi

    fasterq-dump \\
        --threads ${task.cpus} \\
        --temp /tmp \\
        ${sra.name}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}

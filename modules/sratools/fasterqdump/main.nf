// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

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
    tuple val(meta), path(output), emit: reads
    path "versions.yml"          , emit: versions

    script:
    def config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    // Paired-end data extracted by fasterq-dump (--split-3 the default) always creates
    // *_1.fastq *_2.fastq files but sometimes also an additional *.fastq file
    // for unpaired reads which we ignore here.
    output = meta.single_end ? '*.fastq' : '*{1,2}.fastq'
    """
    eval "\$(vdb-config -o n NCBI_SETTINGS | sed 's/[" ]//g')"
    if [[ ! -f "\${NCBI_SETTINGS}" ]]; then
        mkdir -p "\$(dirname "\${NCBI_SETTINGS}")"
        printf '${config}' > "\${NCBI_SETTINGS}"
    fi

    fasterq-dump \\
        ${options.args} \\
        --threads $task.cpus \\
        ${sra.name}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}

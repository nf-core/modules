// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUSCO {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::busco=5.2.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/busco:5.2.2--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)
    path(augustus_config)
    val(lineage)

    output:
    tuple val(meta), path("${meta.id}/run_*/full_table.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}/run_*/short_summary.txt"), emit: txt
    path "versions.yml", emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (lineage) options.args += " --lineage_dataset $lineage"
    """
    # Ensure the input is uncompressed
    gzip -cdf $fasta > __UNCOMPRESSED_FASTA_FILE__
    # Copy the image's AUGUSTUS config directory if it was not provided to the module
    [ ! -e augustus_config ] && cp -a /usr/local/config augustus_config
    AUGUSTUS_CONFIG_PATH=augustus_config \\
    busco \\
        $options.args \\
        --augustus \\
        --cpu $task.cpus \\
        --in __UNCOMPRESSED_FASTA_FILE__ \\
        --out $meta.id

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}

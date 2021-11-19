include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DASTOOL_SCAFFOLDS2BIN {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::das_tool=1.1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/das_tool:1.1.3--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/das_tool:1.1.3--r41hdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)
    val(extension)

    output:
    tuple val(meta), path("*.tsv"), emit: scaffolds2bin
    path "versions.yml"                         , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def file_extension = extension ? extension : "fasta"

    """
    gunzip -f *.${file_extension}.gz

    Fasta_to_Scaffolds2Bin.sh \\
        $options.args \\
        -i . \\
        -e $file_extension \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( DAS_Tool --version 2>&1 | grep "DAS Tool" | sed 's/DAS Tool version //' )
    END_VERSIONS
    """
}

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '3.23'

process MUMMER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::mummer=3.23" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mummer:3.23--pl5262h1b792b2_12"
    } else {
        container "quay.io/biocontainers/mummer:3.23--pl5262h1b792b2_12"
    }

    input:
    tuple val(meta), path(ref), path(query)

    output:
    tuple val(meta), path("*.coords"), emit: coords
    path "versions.yml"              , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def is_compressed_ref = ref.getName().endsWith(".gz") ? true : false
    def fasta_name_ref = ref.getName().replace(".gz", "")

    def is_compressed_query = query.getName().endsWith(".gz") ? true : false
    def fasta_name_query = query.getName().replace(".gz", "")
    """
    if [ "$is_compressed_ref" == "true" ]; then
        gzip -c -d $ref > $fasta_name_ref
    fi
    if [ "$is_compressed_query" == "true" ]; then
        gzip -c -d $query > $fasta_name_query
    fi
    mummer \\
        $options.args \\
        $fasta_name_ref \\
        $fasta_name_query \\
        > ${prefix}.coords

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo $VERSION )
    END_VERSIONS
    """
}

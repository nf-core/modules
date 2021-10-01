// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALLELECOUNTER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::cancerit-allelecount=4.3.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cancerit-allelecount:4.3.0--h41abebc_0"
    } else {
        container "quay.io/biocontainers/cancerit-allelecount:4.3.0--h41abebc_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path loci

    output:
    tuple val(meta), path("*.alleleCount"), emit: allelecount
    path "versions.yml"                   , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    alleleCounter \\
        $options.args \\
        -l $loci \\
        -b $bam \\
        -o ${prefix}.alleleCount

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(alleleCounter --version)
    END_VERSIONS
    """
}

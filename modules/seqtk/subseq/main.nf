// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQTK_SUBSEQ {
    tag '$sequences'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
    } else {
        container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
    }

    input:
    path sequences
    path filter_list

    output:
    path "*.gz"         , emit: sequences
    path "versions.yml" , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ?: ''
    def ext = "fa"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        ext = "fq"
    }
    """
    seqtk \\
        subseq \\
        $options.args \\
        $sequences \\
        $filter_list | \\
        gzip --no-name > ${sequences}${prefix}.${ext}.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(seqtk 2>&1 | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

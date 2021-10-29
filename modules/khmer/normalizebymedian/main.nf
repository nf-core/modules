// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KHMER_NORMALIZEBYMEDIAN {
    tag "${name}"
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::khmer=3.0.0a3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2"
    } else {
        container "quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2"
    }

    input:
    path pe_reads
    path se_reads
    val  name

    output:
    path "${name}.fastq.gz", emit: reads
    path "versions.yml"    , emit: versions

    script:
    pe_args = pe_reads ? "--paired" : ""
    se_args = se_reads ? "--unpaired-reads ${se_reads}" : ""
    files   = pe_reads ? pe_reads : se_reads

    """
    normalize-by-median.py -M ${task.memory.toGiga()}e9 --gzip ${options.args} -o ${name}.fastq.gz ${pe_args} ${se_args} ${files}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}

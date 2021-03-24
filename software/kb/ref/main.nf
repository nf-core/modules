// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KB_REF {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::kb-python==0.25.1--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kb-python:0.25.1--py_0"
    } else {
        container "quay.io/biocontainers/kb-python:0.25.1--py_0"
    }

    input:
    path(fasta)
    path(gtf)

    output:
    path "kb_ref_index" ,   emit: idx
    path "t2g"          ,   emit: t2g
    path "cdna"          ,  emit: cdna
    path "*.version.txt",   emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    kb \\
    ref \\
    $options.args \\
    -i kb_ref_index \\
    -g t2g \\
    -f1 cdna \\
    $fasta \\
    $gtf

    echo \$(kb 2>&1) | sed 's/^kb_python //; s/Usage.*\$//' > ${software}.version.txt
    """
}

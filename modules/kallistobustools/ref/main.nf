// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KALLISTOBUSTOOLS_REF {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::kb-python=0.26.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kb-python:0.26.3--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/kb-python:0.26.3--pyhdfd78af_0"
    }

    input:
    path fasta
    path gtf
    val  workflow

    output:
    path "versions.yml"   , emit: versions
    path "kb_ref_out.idx" , emit: index
    path "t2g.txt"        , emit: t2g
    path "cdna.fa"        , emit: cdna
    path "intron.fa"      , optional:true, emit: intron
    path "cdna_t2c.txt"   , optional:true, emit: cdna_t2c
    path "intron_t2c.txt" , optional:true, emit: intron_t2c

    script:
    if (workflow == "standard") {
        """
        kb \\
            ref \\
            -i kb_ref_out.idx \\
            -g t2g.txt \\
            -f1 cdna.fa \\
            --workflow $workflow \\
            $fasta \\
            $gtf

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
        END_VERSIONS
        """
    } else {
        """
        kb \\
            ref \\
            -i kb_ref_out.idx \\
            -g t2g.txt \\
            -f1 cdna.fa \\
            -f2 intron.fa \\
            -c1 cdna_t2c.txt \\
            -c2 intron_t2c.txt \\
            --workflow $workflow \\
            $fasta \\
            $gtf

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
        END_VERSIONS
        """
    }
}

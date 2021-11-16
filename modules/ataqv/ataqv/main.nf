include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ATAQV_ATAQV {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ataqv=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ataqv:1.2.1--py39ha23c084_2"
    } else {
        container "quay.io/biocontainers/ataqv:1.2.1--py36hfdecbe1_2"
    }

    input:
    tuple val(meta), path(bam), path(bai), path(peak_file)
    val organism
    path tss_file
    path excl_regs_file
    path autosom_ref_file

    output:
    tuple val(meta), path("*.ataqv.json"), emit: json
    tuple val(meta), path("*.problems")  , emit: problems, optional: true
    path "versions.yml"                  , emit: versions

    script:
    def prefix      = options.suffix   ? "${meta.id}${options.suffix}"                  : "${meta.id}"
    def peak        = peak_file        ? "--peak-file $peak_file"                       : ''
    def tss         = tss_file         ? "--tss-file $tss_file"                         : ''
    def excl_regs   = excl_regs_file   ? "--excluded-region-file $excl_regs_file"       : ''
    def autosom_ref = autosom_ref_file ? "--autosomal-reference-file $autosom_ref_file" : ''
    """
    ataqv \\
        $options.args \\
        $peak \\
        $tss \\
        $excl_regs \\
        $autosom_ref \\
        --metrics-file "${prefix}.ataqv.json" \\
        --threads $task.cpus \\
        --name $prefix \\
        $organism \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( ataqv --version )
    END_VERSIONS
    """
}

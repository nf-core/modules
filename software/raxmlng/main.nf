// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RAXMLNG {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::raxml-ng=1.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/raxml-ng:1.0.2--h7447c1b_0"
    } else {
        container "quay.io/biocontainers/raxml-ng:1.0.2--h7447c1b_0"
    }

    input:
    path alignment

    output:
    path "*.raxml.bestTree", emit: phylogeny
    path "*.raxml.support" , optional:true, emit: phylogeny_bootstrapped
    path "*.version.txt"   , emit: version

    script:
    def software = getSoftwareName(task.process)

    if (options.args.contains('--bs-trees')) {
        options.args = "--all ${options.args}"
    }
    """
    raxml-ng \\
        $options.args \\
        --msa $alignment \\
        --threads $task.cpus \\
        --prefix output

    echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//' > ${software}.version.txt
    """
}

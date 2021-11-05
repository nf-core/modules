// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DAMAGEPROFILER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::damageprofiler=1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/damageprofiler:1.1--hdfd78af_2"
    } else {
        container "quay.io/biocontainers/damageprofiler:1.1--hdfd78af_2"
    }

    input:
    tuple val(meta), path(bam)
    path fasta
    path fai
    path specieslist

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "versions.yml"              , emit: versions

    script:
    def software     = getSoftwareName(task.process)
    prefix           = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def reference    = fasta ? "-r $fasta" : ""
    def species_list = specieslist ? "-sf $specieslist" : ""

    """
    damageprofiler \\
    -i $bam \\
    -o $prefix/ \\
    $options.args \\
    $reference \\
    $species_list

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(damageprofiler -v | sed 's/^DamageProfiler v//')
    END_VERSIONS
    """

}

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_GENOMICSDBIMPORT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(vcf), path(tbi), path(intervalfile), val(intervalval), path(wspace)
    val run_intlist
    val run_updatewspace
    val input_map

    output:
    tuple val(meta), path("${prefix}")      , optional:true, emit: genomicsdb
    tuple val(meta), path("$updated_db")    , optional:true, emit: updatedb
    tuple val(meta), path("*.interval_list"), optional:true, emit: intervallist
    path "versions.yml"                                    , emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // settings for running default create gendb mode
    inputs_command = input_map ? "--sample-name-map ${vcf[0]}" : "${'-V ' + vcf.join(' -V ')}"
    dir_command = "--genomicsdb-workspace-path ${prefix}"
    intervals_command = intervalfile ? " -L ${intervalfile} " : " -L ${intervalval} "

    // settings changed for running get intervals list mode if run_intlist is true
    if (run_intlist) {
        inputs_command = ''
        dir_command = "--genomicsdb-update-workspace-path ${wspace}"
        intervals_command = "--output-interval-list-to-file ${prefix}.interval_list"
    }

    // settings changed for running update gendb mode. inputs_command same as default, update_db forces module to emit the updated gendb
    if (run_updatewspace) {
        dir_command = "--genomicsdb-update-workspace-path ${wspace}"
        intervals_command = ''
        updated_db = wspace.toString()
    }

    """
    gatk GenomicsDBImport \\
        $inputs_command \\
        $dir_command \\
        $intervals_command \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

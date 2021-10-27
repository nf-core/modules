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

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(vcf), path(tbi), path(wspace), path(intervalfile), val(intervalval)
    val run_intlist
    val run_updatewspace
    val input_map

    output:
    tuple val(meta), path("*_genomicsdb")   , optional:true, emit: genomicsdb
    tuple val(meta), path("*.interval_list"), optional:true, emit: intervallist
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def inputs_list = []
    def inputs_command = ''
    def dir_command = ''
    def intervals_command = ''

    if(run_intlist){
        inputs_command = ''
        dir_command = "--genomicsdb-update-workspace-path ${wspace}"
        intervals_command = "--output-interval-list-to-file ${prefix}.interval_list"
    } else {
        if(input_map){
            inputs_command = "--sample-name-map ${vcf[0]}"
        } else {
            vcf.each() {a -> inputs_list.add(" -V " + a)}
            inputs_command = inputs_list.join(' ')
        }

        if(run_updatewspace){
            dir_command = "--genomicsdb-update-workspace-path ${wspace}"
            intervals_command = ''
        } else {
            dir_command = "--genomicsdb-workspace-path ${prefix}_genomicsdb"
            intervals_command = intervalfile ? " -L ${intervalfile} " : " -L ${intervalval} "
        }
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

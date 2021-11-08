// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK3_UNIFIEDGENOTYPER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::gatk=3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk:3.5--hdfd78af_11"
    } else {
        container "quay.io/biocontainers/gatk:3.5--hdfd78af_11"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path fai
    path dict

    output:
    path "*.vcf.gz",        emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    gatk3 -T UnifiedGenotyper \\
        -R $ref \\
        -I $bam \\
        -o ${prefix}.vcf \\
        $options.args

    gzip --no-name ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( gatk3 --version )
    END_VERSIONS
    """
}

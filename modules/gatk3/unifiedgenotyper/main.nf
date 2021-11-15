// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK3_UNIFIEDGENOTYPER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk=3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-aa1d7bddaee5eb6c4cbab18f8a072e3ea7ec3969:f963c36fd770e89d267eeaa27cad95c1c3dbe660-0"
    } else {
        container "quay.io/biocontainers/gatk:3.5--hdfd78af_11"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz") , emit: vcf
    path "versions.yml"               , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    gatk3 -T UnifiedGenotyper \\
        -R $ref \\
        -I $bam \\
        -o ${prefix}.vcf \\
        $options.args

    bgzip --no-name ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( gatk3 --version )
    END_VERSIONS
    """
}

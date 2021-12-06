// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'


process GATK3_UNIFIEDGENOTYPER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk=3.5 bioconda::htslib=1.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-aa1d7bddaee5eb6c4cbab18f8a072e3ea7ec3969:f963c36fd770e89d267eeaa27cad95c1c3dbe660-0' :
        'quay.io/biocontainers/mulled-v2-aa1d7bddaee5eb6c4cbab18f8a072e3ea7ec3969:f963c36fd770e89d267eeaa27cad95c1c3dbe660-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz") , emit: vcf
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk3 -T UnifiedGenotyper \\
        -R $ref \\
        -I $bam \\
        -o ${prefix}.vcf \\
        $args

    bgzip \\
        $args2 \\
        -@ ${task.cpus} \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( gatk3 --version )
    END_VERSIONS
    """
}

process GATK4_MERGEVCFS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs)
    path  ref_dict
    val   use_ref_dict

    output:
    tuple val(meta), path('*.vcf.gz'), emit: vcf
    path  "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"

    // Make list of VCFs to merge
    def input = ""
    for (vcf in vcfs) {
        input += " I=${vcf}"
    }
    def ref = use_ref_dict ? "D=${ref_dict}" : ""
    """
    gatk MergeVcfs \\
        $input \\
        O=${prefix}.vcf.gz \\
        $ref \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

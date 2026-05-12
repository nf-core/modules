process OCTOPUSV_CLEAN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/octopusv:0.3.2--pyhdfd78af_0' :
        'quay.io/biocontainers/octopusv:0.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_arg = fasta ? "-g ${fasta}" : ""
    """
    vcf_in="${vcf}"
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} > ${prefix}.input.vcf
        vcf_in="${prefix}.input.vcf"
    fi

    octopusv clean \\
        \${vcf_in} \\
        ${prefix}.vcf.gz \\
        ${fasta_arg} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    echo "" | gzip > ${prefix}.vcf.gz
    """
}

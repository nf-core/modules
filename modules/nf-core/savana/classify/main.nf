process SAVANA_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/savana:1.3.7--pyhdfd78af_0'
        : 'quay.io/biocontainers/savana:1.3.7--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.classified.vcf")              , emit: classified_vcf
    tuple val(meta), path("${prefix}.classified.somatic.vcf")      , emit: somatic_vcf        , optional: true
    tuple val(meta), path("${prefix}.classified.germline.vcf")     , emit: germline_vcf       , optional: true
    tuple val(meta), path("${prefix}.classified.somatic.bedpe")    , emit: somatic_bedpe      , optional: true
    tuple val(meta), path("${prefix}.classified.{strict,lenient}.vcf"), emit: legacy_vcfs     , optional: true
    tuple val("${task.process}"), val('savana'), eval("python -c \"import importlib.metadata as m; print(m.version('savana'))\""), emit: versions_savana, topic: versions


    script:
    def args = (task.ext.args ?: '').trim()
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf_in="${vcf}"
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} > ${prefix}.input.vcf
        vcf_in="${prefix}.input.vcf"
    fi

    savana classify \\
        --vcf \${vcf_in} \\
        --output ${prefix}.classified.vcf \\
        --somatic_output ${prefix}.classified.somatic.vcf \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.classified.vcf
    touch ${prefix}.classified.somatic.vcf
    touch ${prefix}.classified.somatic.bedpe
    touch ${prefix}.classified.strict.vcf
    touch ${prefix}.classified.lenient.vcf
    """
}

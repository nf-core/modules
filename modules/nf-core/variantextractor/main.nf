process VARIANTEXTRACTOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/variant-extractor:5.1.0--pyh106432d_0':
        'biocontainers/variant-extractor:5.1.0--pyh106432d_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions_variantextractor, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix       = task.ext.prefix ?: "${meta.id}"
    pass_only    = task.ext.args?.contains('--pass-only')       ? 'True'  : 'False'
    ensure_pairs = task.ext.args?.contains('--no-ensure-pairs') ? 'False' : 'True'

    """
    echo ${pass_only}
    echo ${ensure_pairs}
    """

    template 'variantextractor.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.extracted.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        variant-extractor: \$(python3 -c 'import variant_extractor; print(variant_extractor.__version__)')
    END_VERSIONS
    """
}

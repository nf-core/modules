process CNVKIT_EXPORT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::cnvkit=0.9.9 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0':
        'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(cns)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: output
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.args.tokenize(" ")[0]
    """
    cnvkit.py export \\
        $args \\
        $cns \\
        -o ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
    END_VERSIONS
    """
}

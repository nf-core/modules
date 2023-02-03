process CNVKIT_CALL {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::cnvkit=0.9.9 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0':
        'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(cns), path(vcf)

    output:
    tuple val(meta), path("*.cns"), emit: cns
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_cmd = vcf ? "-v $vcf" : ""
    """
    cnvkit.py call \\
        $cns \\
        $vcf_cmd \\
        $args \\
        -o ${prefix}.cns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
    END_VERSIONS
    """
}

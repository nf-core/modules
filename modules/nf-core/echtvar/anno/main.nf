process ECHTVAR_ANNO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8fbd0fb4d4e5dc62fbec999fb80befec797d83a0b3b70f39b6baaa73714438df/data':
        'community.wave.seqera.io/library/echtvar_htslib:d37d5e1f4106f9c3' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(annotation_files)
    val(output_suffix)

    output:
    tuple val(meta), path("*.{vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val("${task.process}"), val('echtvar'), eval("echtvar --version | sed 's/echtvar //g'"), topic: versions, emit: versions_echtvar

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = output_suffix ?: 'vcf.gz'
    def annotations = annotation_files ? "-e ${annotation_files.join(' -e ')}" : ''
    """
    echtvar \\
        anno \\
        $args \\
        $annotations \\
        $vcf \\
        ${prefix}.${suffix}

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = output_suffix ?: 'vcf.gz'
    def create_output = suffix.endsWith('.gz') ? "echo '' | bgzip -c > ${prefix}.${suffix}" : "touch ${prefix}.${suffix}"
    """
    echo $args

    ${create_output}
    """
}

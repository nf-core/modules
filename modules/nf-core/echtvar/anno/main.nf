process ECHTVAR_ANNO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87b75cb9e32b89261e8cbdca40b219a5d58fc78ebf92d5ca97c7ca23da1b9517/data':
        'community.wave.seqera.io/library/echtvar:0.2.4--e59eba33636e3aab' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(annotation_files)
    val(output_suffix)

    output:
    tuple val(meta), path("*.{vcf.gz,bcf}"), emit: vcf
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
    """
    echo $args

    echo "" | gzip > ${prefix}.${suffix}
    """
}

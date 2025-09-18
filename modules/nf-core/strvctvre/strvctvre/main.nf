process STRVCTVRE_STRVCTVRE {
    tag "$meta.id"
    label 'process_low'

    // When updating the version here, don't forget to update the hardcoded version in script and stub sections
    // Version in help message does not match the actual version of the tool
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/9584eeb6569a511be29d0a07bf80103d59d38715ddb971dddeca0bc72aec41d3/data':
        'community.wave.seqera.io/library/liftover_strvctvre:5fec172b808cc48e' }"

    input:
    tuple val(meta), path(sv_file), path(sv_file_index), val(assembly)
    tuple val(meta2), path(phylop)
    tuple val(meta3), path(data_directory)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf, optional: true
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.10'
    def format = ''
    if (sv_file.name.endsWith('.vcf') || sv_file.name.endsWith('.vcf.gz')) {
        format = 'vcf'
    } else if (sv_file.name.endsWith('.bed')) {
        format = 'bed'
    } else {
        error("Input structural variants file must be in VCF or BED format")
    }
    if (!['GRCh38', 'GRCh37'].contains(assembly)) {
        error("Assembly must be either 'GRCh37' or 'GRCh38'")
    }
    """
    StrVCTVRE.py \\
        --input ${sv_file} \\
        --format ${format} \\
        --phyloP ${phylop} \\
        --assembly ${assembly} \\
        --liftover liftover_hg19_to_hg38_public.py \\
        --output ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strvctvre: '$VERSION'
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.10'
    def format = ''
    if (sv_file.name.endsWith('.vcf') || sv_file.name.endsWith('.vcf.gz')) {
        format = 'vcf'
    } else if (sv_file.name.endsWith('.bed')) {
        format = 'bed'
    } else {
        error("Input structural variants file must be in VCF or BED format")
    }
    if (!['GRCh38', 'GRCh37'].contains(assembly)) {
        error("Assembly must be either 'GRCh37' or 'GRCh38'")
    }
    """
    touch ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strvctvre: '$VERSION'
    END_VERSIONS
    """
}

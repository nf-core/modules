// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VCF2MAF {

    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version strings below when bumping container versions.
    // ensembl-vep 113.4 still reports "113.0" on the CLI.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/dd6d05067c8171db981b0a2e914561ba7ec4340719ffa8e4d39fcef7004fb230/data'
        : 'community.wave.seqera.io/library/ensembl-vep_vcf2maf:032bb9f304259c9c'}"

    input:
    tuple val(meta), path(vcf) // Use an uncompressed VCF file!
    path fasta                 // Required
    path vep_cache             // Required for VEP running. A default of /.vep is supplied.

    output:
    tuple val(meta), path("*.maf"), emit: maf
    tuple val("${task.process}"), val('vcf2maf'), val("1.6.22"), emit: versions_vcf2maf, topic: versions
    tuple val("${task.process}"), val('ensemblvep'), val("113.4"), emit: versions_ensemblvep, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def vep_cache_cmd = vep_cache       ? "--vep-data $vep_cache" : "" // If VEP is present, it will find it and add it to commands otherwise blank
    """
    if [ "$vep_cache" ]; then
        VEP_CMD="--vep-path \$(dirname \$(type -p vep))"
    else
        VEP_CMD=""
    fi

    vcf2maf.pl \\
        $args \\
        \$VEP_CMD \\
        $vep_cache_cmd \\
        --ref-fasta $fasta \\
        --input-vcf $vcf \\
        --output-maf ${prefix}.maf
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.maf
    """
}

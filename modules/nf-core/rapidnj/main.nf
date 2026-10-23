// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RAPIDNJ {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db61bb5e5b32efe506e3545e098ca416db74be8ef818b1596f741eca952b03b1/data'
        : 'community.wave.seqera.io/library/rapidnj_biopython:4a538e3d630e58e7'}"

    input:
    path alignment

    output:
    path "*.sth"       , emit: stockholm_alignment
    path "*.tre"       , emit: phylogeny
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('rapidnj'), val('2.3.2'), emit: versions_rapidnj, topic: versions
    tuple val("${task.process}"), val('biopython'), eval('python -c "import Bio; print(Bio.__version__)"'), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python \\
        -c 'from Bio import SeqIO; SeqIO.convert("$alignment", "fasta", "alignment.sth", "stockholm")'

    rapidnj \\
        alignment.sth \\
        $args \\
        -i sth \\
        -c $task.cpus \\
        -x rapidnj_phylogeny.tre
    """

    stub:
    """
    touch alignment.sth
    touch rapidnj_phylogeny.tre
    """
}

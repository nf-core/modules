process RPBP_PREPAREGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("${prefix}")                                                       , emit: index
    tuple val(meta), path("${prefix}/${name}.annotated.bed.gz")                              , emit: transcript_bed
    tuple val(meta), path("${prefix}/transcript-index/${name}.orfs-genomic.annotated.bed.gz"), emit: orfs_genomic_bed
    tuple val(meta), path("${prefix}/transcript-index/${name}.orfs-exons.annotated.bed.gz")  , emit: orfs_exons_bed
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    name   = meta.id ?: 'reference'
    prefix = task.ext.prefix ?: 'rpbp_index'
    """
    mkdir -p ${prefix}/transcript-index ${prefix}/star
    grep '>' ${fasta} | sed 's/>//; s/ .*//' > ${prefix}/star/chrName.txt

    python <<'PYTHON'
import argparse
from rpbp.reference_preprocessing.prepare_rpbp_genome import get_orfs

config = {
    "genome_base_path": "${prefix}",
    "genome_name":      "${name}",
    "gtf":              "${gtf}",
    "fasta":            "${fasta}",
    "star_index":       "${prefix}/star",
}

args = argparse.Namespace(
    do_not_call=False,
    overwrite=False,
    num_cpus=${task.cpus},
    log_file="",
    enable_ext_logging=False,
    log_stdout=False,
    no_log_stderr=False,
    logging_level="WARNING",
    file_logging_level="NOTSET",
    stdout_logging_level="NOTSET",
    stderr_logging_level="NOTSET",
)

get_orfs(config["gtf"], args, config, is_annotated=True, is_de_novo=False)
PYTHON

    rm -rf ${prefix}/star
    """

    stub:
    name   = meta.id ?: 'reference'
    prefix = task.ext.prefix ?: 'rpbp_index'
    """
    mkdir -p ${prefix}/transcript-index
    echo "" | gzip > ${prefix}/${name}.annotated.bed.gz
    echo "" | gzip > ${prefix}/transcript-index/${name}.orfs-genomic.annotated.bed.gz
    echo "" | gzip > ${prefix}/transcript-index/${name}.orfs-exons.annotated.bed.gz
    """
}

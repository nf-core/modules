// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ARTIC_MINION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6019bce7e1652ddff248bbfe2a6d9fc1c668137b731cd5ed772e690a7d5ddca/data' :
        'community.wave.seqera.io/library/artic:1.12.0--ac55155577657f92' }"

    input:
    tuple val(meta), path(fastq)
    tuple val(meta2), path(model_dir), val(model)
    tuple val(meta3), path(fasta), path(bed)
    path hdf5_plugin_path

    output:
    tuple val(meta), path("${prefix}.*")                              , emit: results
    tuple val(meta), path("${prefix}.sorted.bam")                     , emit: bam
    tuple val(meta), path("${prefix}.sorted.bam.bai")                 , emit: bai
    tuple val(meta), path("${prefix}.primertrimmed.rg.sorted.bam")    , emit: bam_primertrimmed
    tuple val(meta), path("${prefix}.primertrimmed.rg.sorted.bam.bai"), emit: bai_primertrimmed
    tuple val(meta), path("${prefix}.consensus.fasta")                , emit: fasta
    tuple val(meta), path("${prefix}.pass.vcf.gz")                    , emit: vcf
    tuple val(meta), path("${prefix}.pass.vcf.gz.tbi")                , emit: tbi
    tuple val(meta), path("*.json")                                   , emit: json, optional:true
    tuple val("${task.process}"), val('artic'), eval("artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//'"), topic: versions, emit: versions_artic

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def model_dir_val   = model_dir ?: "\$(which artic | sed 's/artic/models/')"
    def hd5_plugin_path = hdf5_plugin_path?: "/usr/local/lib/python3.6/site-packages/ont_fast5_api/vbz_plugin"
    """
    export HDF5_PLUGIN_PATH=${hd5_plugin_path}

    artic \\
        minion \\
        ${args} \\
        --threads ${task.cpus} \\
        --read-file ${fastq} \\
        --bed ${bed} \\
        --ref ${fasta} \\
        --model-dir ${model_dir_val} \\
        --model ${model} \\
        ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip() { pigz; }

    touch ${prefix}.1.vcf
    touch ${prefix}.2.vcf
    touch ${prefix}.unmatched.vcf

    touch ${prefix}.alignreport.tsv
    touch ${prefix}.amplicon_depths.tsv
    touch ${prefix}.amplicons.tsv

    touch ${prefix}.consensus.fasta
    touch ${prefix}.coverage_mask.txt

    touch ${prefix}.fail.vcf
    touch ${prefix}.ignore.vcf

    touch ${prefix}.merged.vcf
    echo "" | gzip > ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.gz.tbi

    touch ${prefix}.minion.log.txt
    touch ${prefix}.overlaps.tsv

    echo "" | gzip > ${prefix}.normalised.vcf.gz
    touch ${prefix}.normalised.vcf.gz.tbi

    touch ${prefix}.pass.vcf
    echo "" | gzip > ${prefix}.pass.vcf.gz
    touch ${prefix}.pass.vcf.gz.tbi

    touch ${prefix}.preconsensus.fasta
    touch ${prefix}.preconsensus.fasta.fai
    touch ${prefix}.pre-normalisation.coverage.tsv
    touch ${prefix}.post-normalisation.coverage.tsv

    touch ${prefix}.primers.vcf
    touch ${prefix}.primersitereport.txt
    touch ${prefix}.primer_pairs.tsv
    touch ${prefix}.primertrimmed.rg.sorted.bam
    touch ${prefix}.primertrimmed.rg.sorted.bam.bai

    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch ${prefix}.split_reads.tsv
    touch ${prefix}.summary.tsv
    """
}

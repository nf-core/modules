process RASTAIR2_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50c87ecbec7ea76f99ff4ad84e16f083c0b74651bf813d9158b1004234f6feca/data'
        : 'community.wave.seqera.io/library/htslib_rastair:4ffd26b40813e268'}"

    input:
    tuple val(meta), path(bam), path(index), path(mbias_cutoffs)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.bed.gz"), emit: bed
    tuple val(meta), path("*.bed.gz.tbi"), emit: bed_tbi
    tuple val("${task.process}"), val('rastair'), eval("rastair --version | sed 's/^rastair //'"), emit: versions_rastair, topic: versions
    tuple val("${task.process}"), val('tabix'), eval("tabix --version 2>&1 | sed -n '1s/^tabix (htslib) //p'"), emit: versions_tabix, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Single-end reads carry no pairing flags: strand comes from the alignment direction alone.
    def unpaired = meta.single_end ? '--unpaired' : ''
    // Bases to ignore at read ends, as suggested by rastair2/mbias: line 1 is --nOT, line 2 --nOB
    def cutoffs_files = mbias_cutoffs instanceof List ? mbias_cutoffs : [mbias_cutoffs]
    if (cutoffs_files.size() > 1) {
        error("rastair2/call takes one M-bias cutoffs file per sample, got ${cutoffs_files.size()}: restrict rastair2/mbias to one contig with --region")
    }
    def trim = cutoffs_files ? "--nOT \$(sed -n 1p ${mbias_cutoffs}) --nOB \$(sed -n 2p ${mbias_cutoffs})" : ''
    """
    rastair call \\
        ${args} \\
        ${unpaired} \\
        ${trim} \\
        --threads ${task.cpus} \\
        --fasta-file ${fasta} \\
        --vcf ${prefix}.vcf.gz \\
        --bed ${prefix}.bed.gz \\
        ${bam}

    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi
    """
}

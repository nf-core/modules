process RASTAIR2_MBIAS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50c87ecbec7ea76f99ff4ad84e16f083c0b74651bf813d9158b1004234f6feca/data'
        : 'community.wave.seqera.io/library/htslib_rastair:4ffd26b40813e268'}"

    input:
    tuple val(meta), path(bam), path(index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.mbias.html"), emit: html
    tuple val(meta), path("*.cutoffs.txt"), emit: cutoffs, optional: true
    tuple val("${task.process}"), val('rastair'), eval("rastair --version | sed 's/^rastair //'"), emit: versions_rastair, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The report tells read 1 from read 2 and OT from OB with read-pair flags, so it needs paired-end reads.
    if (meta.single_end) {
        error("rastair2 mbias needs paired-end reads: it assigns strands and read ends from read-pair flags")
    }
    // `rastair mbias --bam` resolves symlinks and misses the staged .bai/.fai, so per-read runs here
    """
    rastair per-read \\
        ${args} \\
        --threads ${task.cpus} \\
        --fasta-file ${fasta} \\
        --bed ${prefix}.per_read.bed.gz \\
        ${bam}

    # rastair looks for its R scripts next to the binary, but conda installs them under share/
    rastair mbias \\
        ${args2} \\
        --r-script-dir \$(dirname \$(which rastair))/../share/rastair/scripts \\
        --no-vbias \\
        --no-gc \\
        --threads ${task.cpus} \\
        --bed ${prefix}.per_read.bed.gz \\
        --output-prefix .

    rm ${prefix}.per_read.bed.gz ${prefix}.per_read.bed.gz.tbi
    mv qc_report.html ${prefix}.mbias.html
    for cutoffs in *_cutoffs.txt; do
        [ -e "\$cutoffs" ] || continue
        mv "\$cutoffs" "${prefix}.\${cutoffs%_cutoffs.txt}.cutoffs.txt"
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mbias.html
    touch ${prefix}.chr1.cutoffs.txt
    """
}

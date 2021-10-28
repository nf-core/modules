process BAMCMP {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bamcmp" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://crukmi/bamcmp:2.0.0"
    } else {
        container "quay.io/biocontainers/bamcmp"
    }

    input:
    tuple val(meta), path(bam)

    output:
    path "*.bam", emit: bam
    path "versions.yml"          , emit: versions

    script:

    """
    bamcmp -s "as" ${bam} \\
    -1 ${bam[0]} \\
    -2 ${bam[1]} \\
    -A ${meta.id}_contaminationFiltered.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}

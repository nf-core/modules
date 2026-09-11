process GATK4_CONDENSEDEPTHEVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01'}"

    input:
    tuple val(meta), path(depth_evidence), path(depth_evidence_index)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.rd.txt.gz"), emit: condensed_evidence
    tuple val(meta), path("*.rd.txt.gz.tbi"), emit: condensed_evidence_index
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (depth_evidence == "${prefix}.rd.txt.gz") {
        error("File name collision - Please specify a different prefix.")
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CondenseDepthEvidence] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CondenseDepthEvidence \\
        --depth-evidence ${depth_evidence} \\
        --output ${prefix}.rd.txt.gz \\
        --reference ${fasta} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.rd.txt.gz
    touch ${prefix}.rd.txt.gz.tbi
    """
}

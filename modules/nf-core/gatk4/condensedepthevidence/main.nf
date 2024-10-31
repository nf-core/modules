process GATK4_CONDENSEDEPTHEVIDENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0':
        'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0' }"

    input:
    tuple val(meta), path(depth_evidence), path(depth_evidence_index)
    path(fasta)
    path(fasta_fai)
    path(dict)

    output:
    tuple val(meta), path("*.rd.txt.gz")    , emit: condensed_evidence
    tuple val(meta), path("*.rd.txt.gz.tbi"), emit: condensed_evidence_index
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (depth_evidence == "${prefix}.rd.txt.gz"){
        error("File name collision - Please specify a different prefix.")
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK CondenseDepthEvidence] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CondenseDepthEvidence \\
        --depth-evidence ${depth_evidence} \\
        --output ${prefix}.rd.txt.gz \\
        --reference ${fasta} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.rd.txt.gz
    touch ${prefix}.rd.txt.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

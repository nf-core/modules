process GATK4_SVCLUSTER {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs), path(indices)
    path ploidy_table
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: clustered_vcf
    tuple val(meta), path("*.vcf.gz.tbi")   , emit: clustered_vcf_index
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input = vcfs.collect({"--variant ${it}"}).join(" ")

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK SVCluster] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SVCluster \\
        --output ${prefix}.vcf.gz \\
        --ploidy-table ${ploidy_table} \\
        ${input} \\
        --reference ${fasta} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

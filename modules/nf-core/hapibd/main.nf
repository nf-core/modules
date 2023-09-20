process HAPIBD {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::hap-ibd=1.0.rev20May22.818"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap-ibd:1.0.rev20May22.818--hdfd78af_0':
        'biocontainers/hap-ibd:1.0.rev20May22.818--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path(map)
    path(exclude)


    output:
    tuple val(meta), path("*.hbd.gz"), emit: hbd
    tuple val(meta), path("*.ibd.gz"), emit: ibd
    tuple val(meta), path("*.log")   , emit: log
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def excludesamples_command = exclude ? "excludesamples=$exclude" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[hapibd] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    hap-ibd -Xmx${avail_mem}M \\
        gt=${vcf} \\
        map=${map} \\
        out=${prefix} \\
        $args \\
        ${excludesamples_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hapibd: \$(hap-ibd 2>&1 |head -n1 | sed 's/^hap-ibd.jar  \\[ version //; s/, /rev/; s/ \\]//')
    END_VERSIONS
    """
}

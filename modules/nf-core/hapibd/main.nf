process HAPIBD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap-ibd:1.0.rev20May22.818--hdfd78af_0':
        'quay.io/biocontainers/hap-ibd:1.0.rev20May22.818--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path(map)
    path(exclude)


    output:
    tuple val(meta), path("*.hbd.gz"), emit: hbd
    tuple val(meta), path("*.ibd.gz"), emit: ibd
    tuple val(meta), path("*.log")   , emit: log
    tuple val("${task.process}"), val('hapibd'), eval("hap-ibd 2>&1 | sed '1!d;s/^.* version //;s/,.*//'"), topic: versions, emit: versions_hapibd


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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[hapibd] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    touch ${prefix}.log
    echo "" | gzip > ${prefix}.hbd.gz
    echo "" | gzip > ${prefix}.ibd.gz
    """
}

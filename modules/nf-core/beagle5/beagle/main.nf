process BEAGLE5_BEAGLE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/beagle:5.5_27Feb25.75f--hdfd78af_0'
        : 'biocontainers/beagle:5.5_27Feb25.75f--hdfd78af_0'}"

    input:
    // Including `val(region)` to prevent errors with multi-chromosome VCFs and single-chromosome reference panels.
    // This enhances clarity and simplifies implementation in the subworkflow.
    tuple val(meta), path(vcf), path(vcf_index), path(refpanel), path(refpanel_index), path(genmap), path(exclsamples), path(exclmarkers), val(region)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.log")   , emit: log
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.bglout"

    def ref_command = refpanel ? "ref=${refpanel}" : ""
    def map_command = genmap   ? "map=${genmap}"   : ""
    def region_cmd  = region   ? "chrom=${region}" : ""

    def excludesamples_command = exclsamples ? "excludesamples=${exclsamples}" : ""
    def excludemarkers_command = exclmarkers ? "excludemarkers=${exclmarkers}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[beagle] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    beagle -Xmx${avail_mem}M \\
        gt=${vcf} \\
        out=${prefix} \\
        ${args} \\
        ${ref_command} \\
        ${map_command} \\
        ${region_cmd} \\
        ${excludesamples_command} \\
        ${excludemarkers_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: \$(beagle 2>&1 |head -n1 | sed -rn 's/beagle\\.(.*)\\.jar \\(version (.*)\\)/\\2rev\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.bglout"
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: \$(beagle 2>&1 |head -n1 | sed -rn 's/beagle\\.(.*)\\.jar \\(version (.*)\\)/\\2rev\\1/p')
    END_VERSIONS
    """
}

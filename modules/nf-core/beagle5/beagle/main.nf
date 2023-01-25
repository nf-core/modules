process BEAGLE5_BEAGLE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::beagle=5.2_21Apr21.304"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beagle:5.2_21Apr21.304--hdfd78af_0':
        'quay.io/biocontainers/beagle:5.2_21Apr21.304--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path(refpanel)
    path(genmap)
    path(exclsamples)
    path(exclmarkers)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.log")        , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.bglout"
    def ref_command = refpanel ? "ref=$refpanel" : ""
    def map_command = genmap ? "map=$genmap" : ""
    def excludesamples_command = exclsamples ? "excludesamples=$exclsamples" : ""
    def excludemarkers_command = exclmarkers ? "excludemarkers=$exclmarkers" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[beagle] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    beagle -Xmx${avail_mem}g \\
        gt=${vcf} \\
        out=${prefix} \\
        $args \\
        ${ref_command} \\
        ${map_command} \\
        ${excludesamples_command} \\
        ${excludemarkers_command} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: \$(beagle 2>&1 |head -n1 | sed -rn 's/beagle\\.(.*)\\.jar \\(version (.*)\\)/\\2rev\\1/p')
    END_VERSIONS
    """
}

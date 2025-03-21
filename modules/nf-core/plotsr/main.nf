process PLOTSR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plotsr:1.1.1--pyh7cba7a3_0':
        'biocontainers/plotsr:1.1.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(syri)
    path(fastas)
    val(genomes)
    path(bedpe)
    path(markers)
    path(tracks)
    path(chrord)
    path(chrname)

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def syri_arg    = ( syri instanceof List )  ? syri.collect { "--sr $it" }.join(' ') : "--sr $syri"
    def genomes_txt = genomes.tokenize('\n').withIndex().collect { line, index -> index == 0 ? line : "${fastas[index-1]}\t${line}" }.join('\\n')
    def bedpe_arg   = bedpe     ? "--bedpe $bedpe"      : ''
    def markers_arg = markers   ? "--markers $markers"  : ''
    def tracks_arg  = tracks    ? "--tracks $tracks"    : ''
    def chrord_arg  = chrord    ? "--chrord $chrord"    : ''
    def chrname_arg = chrname   ? "--chrname $chrname"  : ''
    """
    echo -e "${genomes_txt}" \\
        > ${prefix}.genomes.txt

    plotsr \\
        $syri_arg \\
        --genomes ${prefix}.genomes.txt \\
        $bedpe_arg \\
        $markers_arg \\
        $tracks_arg \\
        $chrord_arg \\
        $chrname_arg \\
        $args \\
        -o ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotsr: \$(plotsr --version)
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def syri_arg    = ( syri instanceof List )  ? syri.collect { "--sr $it" }.join(' ') : "--sr $syri"
    def genomes_txt = genomes.tokenize('\n').withIndex().collect { line, index -> index == 0 ? line : "${fastas[index-1]}\t${line}" }.join('\\n')
    def bedpe_arg   = bedpe     ? "--bedpe $bedpe"      : ''
    def markers_arg = markers   ? "--markers $markers"  : ''
    def tracks_arg  = tracks    ? "--tracks $tracks"    : ''
    def chrord_arg  = chrord    ? "--chrord $chrord"    : ''
    def chrname_arg = chrname   ? "--chrname $chrname"  : ''
    """
    echo -e "${genomes_txt}" \\
        > ${prefix}.genomes.txt

    echo \\
    "plotsr \\
        $syri_arg \\
        --genomes ${prefix}.genomes.txt \\
        $bedpe_arg \\
        $markers_arg \\
        $tracks_arg \\
        $chrord_arg \\
        $chrname_arg \\
        $args \\
        -o ${prefix}.png"

    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotsr: \$(plotsr --version)
    END_VERSIONS
    """
}

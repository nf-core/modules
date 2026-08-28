process STITCHR_THIMBLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b08ffe9ec48ad8e89bab8bfe16127c8de77f9bc973b567e01c61276da25ad6e/data':
        'community.wave.seqera.io/library/pip_stitchr:9b1e4db63c6ec900' }"

    input:
    tuple val(meta), path(samplesheet)
    tuple val(meta2), path(alt_codon_usage)
    tuple val(meta3), path(allele_preferences)
    tuple val(meta4), val(species), path(stitchr_data, stageAs: 'Data')

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: thimble
    tuple val("${task.process}"), val('stitchr'), eval("stitchr --version"), topic: versions, emit: versions_stitchr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // optional files:
    def codon_usage = alt_codon_usage ? "-cu ${alt_codon_usage}" : ""
    def allele_preference = allele_preferences ? "-p ${allele_preferences}" : ""
    """
    # resolve the "Data" package import to the IMGT reference staged by stitchrdl, otherwise not readable (root issues)
    export PYTHONPATH="\$PWD\${PYTHONPATH:+:\$PYTHONPATH}"

    thimble \\
        ${args} \\
        -in ${samplesheet} \\
        -o ${prefix}.tsv \\
        ${codon_usage} \\
        ${allele_preference} \\
        -s ${species.toUpperCase()}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}

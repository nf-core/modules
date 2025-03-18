process MERQURYFK_KATCOMP {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/workspace/modules/.nf-test/tests/310dab7bb6d95f8568312500f359752e/work/df/c123d9edb1a8569c2f3d4f7aaf1b42' :
        'community.wave.seqera.io/library/fastk_merquryfk:ea801837b4afd24b' }"

    input:
    tuple val(meta), path(fastk1_hist), path(fastk1_ktab), path(fastk2_hist), path(fastk2_ktab)

    output:
    tuple val(meta), path("*.fi.{png,pdf}"), emit: filled , optional: true
    tuple val(meta), path("*.ln.{png,pdf}"), emit: line   , optional: true
    tuple val(meta), path("*.st.{png,pdf}"), emit: stacked, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_fk1 = fastk1_ktab.find{ it.toString().endsWith(".ktab") }.getBaseName()
    def input_fk2 = fastk2_ktab.find{ it.toString().endsWith(".ktab") }.getBaseName()
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def FASTK_VERSION   = '1.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.1'
    """
    KatComp \\
        $args \\
        -T$task.cpus \\
        ${input_fk1} \\
        ${input_fk2} \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}

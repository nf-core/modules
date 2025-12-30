process GENMOD_MODELS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/31b331bee43c7ff070bdde5460a4102ba31c3bfb0ee0d70197001ff011036555/data' :
        'community.wave.seqera.io/library/genmod_python:31b2fba4d3b7ba6f' }"

    input:
    tuple val(meta), path(input_vcf), path (fam)
    path (reduced_penetrance)

    output:
    tuple val(meta), path("*_models.vcf"), emit: vcf
    tuple val("${task.process}"), val('genmod'), eval("genmod --version | sed 's/^.*genmod version: //'"), emit: versions_genmod, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def family_file =  fam ? "--family_file ${fam}" : ""
    def pen_file    = reduced_penetrance ? "--reduced_penetrance ${reduced_penetrance}" : ""
    """
    genmod \\
        models \\
        $args \\
        $pen_file \\
        $family_file \\
        --processes ${task.cpus} \\
        --outfile ${prefix}_models.vcf \\
        $input_vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_models.vcf
    """
}

process VCFEXPRESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/056c133c8f4cdb8f8025830e951c7f8b02dbf4c78425cb3a3c1b3bf9e3840ae4/data' :
        'community.wave.seqera.io/library/vcfexpress:0.3.4--bbad2b1ffb0f6492'}"

    input:
    tuple val(meta), path(vcf)
    path(prelude) // optional : empty channel [] if not needed

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val("${task.process}"), val('vcfexpress'), eval("vcfexpress --version | sed 's/vcfexpress //'"), topic: versions, emit: versions_vcfexpress

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_express"
    def suffix = task.ext.suffix ?: 'vcf.gz'
    def lua_prelude = prelude ? "--lua-prelude $prelude" : ''

    """
    vcfexpress filter \
    $args \
    $lua_prelude \
    ${vcf} \
    --output ${prefix}.${suffix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_express"
    def suffix = task.ext.suffix ?: 'vcf.gz'

    """
    echo $args

    echo "" | gzip > ${prefix}.${suffix}
    """
}

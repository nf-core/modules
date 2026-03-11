// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process VCFEXPRESS {
    tag "$meta.id"
    label 'single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/056c133c8f4cdb8f8025830e951c7f8b02dbf4c78425cb3a3c1b3bf9e3840ae4/data' :
        'community.wave.seqera.io/library/vcfexpress:0.3.4--bbad2b1ffb0f6492'}"

    input:
    tuple val(meta), path(vcf)
    path(prelude) // optional : empty channel [] if not needed

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val("${task.process}"), val('vcfexpress'), eval("vcfexpress --version"), topic: versions, emit: versions_vcfexpress

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_express"

    def lua_prelude = prelude ? "--lua-prelude $prelude" : ''

    """
    vcfexpress filter \
    $args \
    $lua_prelude \
    ${vcf} \
    --output ${prefix}.vcf.gz    
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_express"
    
    """
    touch ${prefix}.vcf.gz

    echo $args
    """
}

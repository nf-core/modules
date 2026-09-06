process TIDDIT_SV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3ebf66353bdc536851f786d7427e745c0c02f3951e08fc7828d6b86fefda64ce/data' :
        'community.wave.seqera.io/library/tiddit:3.9.7--33e2b6c6e2f37861' }"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(bwa_index)

    output:
    tuple val(meta), path("${prefix}.vcf")         , emit: vcf
    tuple val(meta), path("${prefix}.ploidies.tab"), emit: ploidy
    tuple val("${task.process}"), val('tiddit'), eval("tiddit | sed -n 's/^usage: tiddit-//; s/ .*//p'"), topic: versions, emit: versions_tiddit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bwa_command = bwa_index ? "[[ -d ${bwa_index} ]] && for i in ${bwa_index}/*; do [[ -f ${fasta} && ! \"\$i\" =~ .*\"${fasta}.\".* ]] && ln -s \$i ${fasta}.\${i##*.} || ln -s \$i .; done" : ""

    """
    $bwa_command

    tiddit \\
        --sv \\
        ${args} \\
        --threads $task.cpus \\
        --bam ${input} \\
        --ref ${fasta} \\
        -o ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" > ${prefix}.vcf
    echo "" > ${prefix}.ploidies.tab
    """
}

process JVARKIT_VCFPOLYX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta),  path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.${extension}"), emit: vcf
    tuple val(meta), path("*.tbi")         , emit: tbi, optional: true
    tuple val(meta), path("*.csi")         , emit: csi, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1  = task.ext.args1 ?: ''
    def args2  = meta.vcfpolyx_args ?: (task.ext.args2 ?: ' --tag POLYX --max-repeats 10 ')
    def args3  = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    extension =     args3.contains("--output-type b") || args3.contains("-Ob") ? "bcf.gz" :
                    args3.contains("--output-type u") || args3.contains("-Ou") ? "bcf" :
                    args3.contains("--output-type z") || args3.contains("-Oz") ? "vcf.gz" :
                    args3.contains("--output-type v") || args3.contains("-Ov") ? "vcf" :
                    "vcf"

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p TMP

    bcftools view -O v ${args1} "${vcf}" |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcfpolyx --reference "${fasta}" ${args2} |\\
        bcftools view --output "${prefix}.${extension}" ${args3}

    rm -rf TMP

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        jvarkit: \$(jvarkit -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.${extension}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        jvarkit: \$(jvarkit -v)
    END_VERSIONS
    """
}

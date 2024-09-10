process JVARKIT_VCFFILTERJDK {
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
    tuple val(meta5), path(code)
    output:
    tuple val(meta), path("*.${extension}"), emit: vcf
    tuple val(meta), path("*.tbi")         , emit: tbi, optional: true
    tuple val(meta), path("*.csi")         , emit: csi, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1       = task.ext.args1 ?: ''
    def args2       = meta.vcffilterjdk_args ?: (task.ext.args2 ?: '')
    def args3       = task.ext.args3 ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def script_file = code?"--script \"${code}\"":""

    extension =     getVcfExtension(args3); /* custom function, see below */

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p TMP

    bcftools view -O v ${args1} "${vcf}" |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcffilterjdk ${script_file}  ${args2} |\\
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


// Custom Function to get VCF extension
String getVcfExtension(String args) {
    return args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
        args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
        args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
        args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
        "vcf";
}

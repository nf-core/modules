process JVARKIT_VCF2TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(regions_file)
    tuple val(meta2), path(pedigree)
    output:
    tuple val(meta), path("*.${extension}"), emit: output
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1        = task.ext.args1 ?: ''
    def args2        = task.ext.args2 ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def ped          = pedigree?"--pedigree \"${pedigree}\"":""
    def regions_file = regions_file? (tbi ? " --regions-file" : " --targets-file")+" \"${regions_file}\" ":""
    extension =     getFileExtension(args2); /* custom function, see below */

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p TMP

    bcftools view ${regions_file} -O v ${args1} "${vcf}" |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcf2table ${ped}  ${args2} > "${prefix}.${extension}"

    rm -rf TMP

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        jvarkit: \$(jvarkit -v)
    END_VERSIONS
    """

    stub:
    def args2        = task.ext.args2 ?: ''
    extension  = getFileExtension(args2); /* custom function, see below */
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
String getFileExtension(String args) {
    return args.contains("--format html") ? "html" : "txt"
}

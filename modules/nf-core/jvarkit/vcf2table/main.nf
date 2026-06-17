process JVARKIT_VCF2TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'quay.io/biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(regions_file)
    tuple val(meta2), path(pedigree)

    output:
    tuple val(meta), path("*.${extension}"), emit: output
    tuple val("${task.process}"), val('jvarkit'), eval("jvarkit -v"), emit: versions_jvarkit, topic: versions
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version |& sed '1!d;s/bcftools //'"), emit: versions_bcftools, topic: versions


    script:
    def args1       = task.ext.args  ?: ''
    def args2       = task.ext.args2 ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def ped         = pedigree?"--pedigree \"${pedigree}\"":""
    def regions_opt = regions_file? (tbi ? " --regions-file" : " --targets-file")+" \"${regions_file}\" ":""
    extension       = getFileExtension(args2); /* custom function, see below */

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p TMP

    bcftools view \\
        ${regions_opt} \\
        -Ov ${args1} \\
        "${vcf}" | \\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcf2table \\
            ${ped} \\
            ${args2} \\
        > "${prefix}.${extension}"

    rm -rf TMP
    """

    stub:
    def args2        = task.ext.args2 ?: ''
    extension  = getFileExtension(args2); /* custom function, see below */
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.${extension}"
    """
}

// Custom Function to get VCF extension
String getFileExtension(String args) {
    return args.contains("--format html") ? "html" : "txt"
}

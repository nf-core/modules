process JVARKIT_VCFPOLYX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'quay.io/biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta),  path(vcf), path(tbi), path(regions_file)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.${extension}"), emit: vcf
    tuple val(meta), path("*.tbi")         , emit: tbi, optional: true
    tuple val(meta), path("*.csi")         , emit: csi, optional: true
    tuple val("${task.process}"), val('jvarkit'), eval("jvarkit -v"), emit: versions_jvarkit, topic: versions
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version |& sed '1!d;s/bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1         = task.ext.args ?: ''
    def args2         = task.ext.args2 ?: ''
    def args3         = task.ext.args3 ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def regions_cmd   = regions_file ? (tbi ? " --regions-file" : " --targets-file") + " '${regions_file}' " : ""

    extension  = getVcfExtension(args3); /* custom function, see below */

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p TMP

    bcftools view -O v \\
        ${regions_cmd} \\
        ${args1} \\
        "${vcf}" |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcfpolyx \\
            --reference "${fasta}" \\
            ${args2} |\\
        bcftools view \\
            --output "${prefix}.${extension}" \\
            ${args3}

    rm -rf TMP
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args3  = task.ext.args3 ?: ''
    extension  = getVcfExtension(args3); /* custom function, see below */
    """
    touch "${prefix}.${extension}"
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

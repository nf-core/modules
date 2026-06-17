process GRAPHTYPER_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphtyper:2.7.2--h7d7f7ad_0':
        'quay.io/biocontainers/graphtyper:2.7.2--h7d7f7ad_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(ref_fai)
    path region_file  // can be empty if --region is supplied to task.ext.args

    output:
    tuple val(meta), path("results/*/*.vcf.gz")    , emit: vcf
    tuple val(meta), path("results/*/*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('graphtyper'), eval("graphtyper --help | tail -1 | sed 's/.* //'"), emit: versions_graphtyper, topic: versions


    script:
    def args          = task.ext.args ?: ''
    def bam_path_text = bam.sort().join('\\n')
    def region_text   = region_file ? "--region_file ${region_file}" : ""
    if (!region_file && ! args.contains("region")) {
        error "GRAPHTYPER_GENOTYPE requires either a region file or a region specified using '--region' in ext.args"
    }
    """
    printf "$bam_path_text" > bam_list.txt
    graphtyper \\
        genotype \\
        ${ref} \\
        ${args} \\
        --sams bam_list.txt \\
        --threads ${task.cpus} \\
        ${region_text}
    """

    stub:
    """
    mkdir -p results/test
    echo "" | gzip > results/test/region1.vcf.gz
    echo "" | gzip > results/test/region2.vcf.gz
    touch results/test/region1.vcf.gz.tbi
    touch results/test/region2.vcf.gz.tbi
    """

}

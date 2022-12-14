process RTGTOOLS_VCFEVAL {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rtg-tools=3.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0':
        'quay.io/biocontainers/rtg-tools:3.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(query_vcf), path(query_vcf_tbi)
    tuple path(truth_vcf), path(truth_vcf_tbi)
    path(truth_regions)
    path(evaluation_regions)
    path(sdf)

    output:
    tuple val(meta), path("**results/{done,progress,*.log}") , emit: logs
    tuple val(meta), path("**tp.vcf.gz"), path("**tp.vcf.gz.tbi") , emit: tp
    tuple val(meta), path("**fn.vcf.gz"), path("**fn.vcf.gz.tbi") , emit: fn
    tuple val(meta), path("**fp.vcf.gz"), path("**fp.vcf.gz.tbi") , emit: fp
    tuple val(meta), path("**baseline.vcf.gz"), path("**baseline.vcf.gz.tbi") , emit: baseline
    tuple val(meta), path("**.tsv.gz")                    , emit: roc
    tuple val(meta), path("**results/summary.txt")        , emit: summary
    tuple val(meta), path("**results/phasing.txt")        , emit: phasing
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_regions = truth_regions ? "--bed-regions=$truth_regions" : ""
    def eval_regions = evaluation_regions ? "--evaluation-regions=$evaluation_regions" : ""
    def truth_index = truth_vcf_tbi ? "" : "rtg index $truth_vcf"
    def query_index = query_vcf_tbi ? "" : "rtg index $query_vcf"
    def avail_mem = task.memory.toGiga() + "G"

    """
    $truth_index
    $query_index

    rtg RTG_MEM=$avail_mem vcfeval \\
        $args \\
        --baseline=$truth_vcf \\
        $bed_regions \\
        $eval_regions \\
        --calls=$query_vcf \\
        --output=${prefix}_results \\
        --template=$sdf \\
        --threads=$task.cpus \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}

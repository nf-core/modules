process RTGTOOLS_VCFEVAL {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::rtg-tools=3.12.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0':
        'quay.io/biocontainers/rtg-tools:3.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(truth_vcf), path(truth_vcf_tbi), path(query_vcf), path(query_vcf_tbi), path(bed)
    path(sdf)

    output:
    tuple val(meta), path("*.txt"), emit: results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--bed-regions=$bed" : ""
    def truth_index = truth_vcf_tbi ? "" : "rtg index $truth_vcf"
    def query_index = query_vcf_tbi ? "" : "rtg index $query_vcf"

    sdf_basename = sdf.getBaseName().replace(".tar","")
    tar_decomp = ""
    if((sdf =~ /.tar.gz\b/).find() == true) {
        tar_decomp = "tar -xzf $sdf"
    }

    """
    $tar_decomp

    $truth_index
    $query_index

    rtg vcfeval \\
        $args \\
        --baseline=$truth_vcf \\
        $regions \\
        --calls=$query_vcf \\
        --output=$prefix \\
        --template=$sdf_basename \\
        --threads=$task.cpus \\
        > ${prefix}_results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}

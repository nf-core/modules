
process SVANALYZER_SVBENCHMARK {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/bioconda/bioconda-recipes/issues/37646
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svanalyzer:0.36--pl526_0':
        'biocontainers/svanalyzer:0.36--pl526_0' }"

    input:
    tuple val(meta), path(test), path(test_tbi), path(truth), path(truth_tbi), path(bed)
    tuple val(meta2),path(fasta)
    tuple val(meta3),path(fai)

    output:
    tuple val(meta), path("*.falsenegatives.vcf.gz"), emit: fns
    tuple val(meta), path("*.falsepositives.vcf.gz"), emit: fps
    tuple val(meta), path("*.distances")            , emit: distances
    tuple val(meta), path("*.log")                  , emit: log
    tuple val(meta), path("*.report")               , emit: report
    tuple val("${task.process}"), val('svanalyzer'), val('0.36'), topic: versions, emit: versions_svanalyzer

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SVANALYZER_SVBENCHMARK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_opt = bed ? "-includebed $bed" : ""

    """
    svanalyzer \\
        benchmark \\
        $args \\
        --ref $fasta \\
        --test $test \\
        --truth $truth \\
        --prefix $prefix \\
        $bed_opt

    bgzip ${args2} --threads ${task.cpus} -c ${prefix}.falsenegatives.vcf > ${prefix}.falsenegatives.vcf.gz
    bgzip ${args2} --threads ${task.cpus} -c ${prefix}.falsepositives.vcf > ${prefix}.falsepositives.vcf.gz
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SVANALYZER_SVBENCHMARK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.falsenegatives.vcf.gz
    echo "" | gzip > ${prefix}.falsepositives.vcf.gz
    touch ${prefix}.distances
    touch ${prefix}.log
    touch ${prefix}.report
    """
}

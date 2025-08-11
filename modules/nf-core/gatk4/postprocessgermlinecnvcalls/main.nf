process GATK4_POSTPROCESSGERMLINECNVCALLS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(calls), path(model), path(ploidy)

    output:
    tuple val(meta), path("*_genotyped_intervals.vcf.gz") , emit: intervals, optional: true
    tuple val(meta), path("*_genotyped_segments.vcf.gz")  , emit: segments, optional: true
    tuple val(meta), path("*_denoised.vcf.gz")            , emit: denoised, optional: true
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def calls_command  = calls   ? calls.collect{"--calls-shard-path $it"}.join(' ')  : ""
    def model_command  = model   ? model.collect{"--model-shard-path $it"}.join(' ')  : ""
    def ploidy_command = ploidy  ? "--contig-ploidy-calls ${ploidy}"                  : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK PostProcessGermlineCnvCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"
    export PYTENSOR_FLAGS="base_compiledir=\$PWD"

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PostprocessGermlineCNVCalls \\
        $calls_command \\
        $model_command \\
        $ploidy_command \\
        $args \\
        --output-genotyped-intervals ${prefix}_genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${prefix}_genotyped_segments.vcf.gz \\
        --output-denoised-copy-ratios ${prefix}_denoised.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_genotyped_intervals.vcf.gz
    echo "" | gzip > ${prefix}_genotyped_segments.vcf.gz
    echo "" | gzip > ${prefix}_denoised.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process GATK4_GERMLINECNVCALLER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(tsv), path(intervals), path(ploidy), path(model)

    output:
    tuple val(meta), path("*-cnv-model/*-calls"), emit: cohortcalls, optional: true
    tuple val(meta), path("*-cnv-model/*-model"), emit: cohortmodel, optional: true
    tuple val(meta), path("*-cnv-calls/*-calls"), emit: casecalls  , optional: true
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals ${intervals}"         : ""
    def ploidy_command    = ploidy    ? "--contig-ploidy-calls ${ploidy}"  : ""
    def model_command     = model     ? "--model ${model}"                 : ""
    def input_list        = tsv.collect{"--input $it"}.join(' ')
    def output_command    = model     ? "--output ${prefix}-cnv-calls"     : "--output ${prefix}-cnv-model"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"
    export PYTENSOR_FLAGS="base_compiledir=\$PWD"
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GermlineCNVCaller \\
        $input_list \\
        $ploidy_command \\
        $output_command \\
        --output-prefix $prefix \\
        $args \\
        $intervals_command \\
        $model_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}-cnv-calls/${prefix}-calls
    mkdir -p ${prefix}-cnv-model/${prefix}-model
    mkdir -p ${prefix}-cnv-model/${prefix}-calls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

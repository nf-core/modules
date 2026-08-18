process GATK4_DETERMINEGERMLINECONTIGPLOIDY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1330cb37b3acde07db0a04befa76973fea72d1d10a79542132e3d77950225af/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:b945148cc5e2a616'}"

    input:
    tuple val(meta), path(counts), path(bed), path(exclude_beds)
    tuple val(meta2), path(ploidy_model)
    path contig_ploidy_table

    output:
    tuple val(meta), path("${prefix}-calls"), emit: calls
    tuple val(meta), path("${prefix}-model"), emit: model, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = bed ? "--intervals ${bed}" : ""
    def exclude = exclude_beds ? exclude_beds.collect { bed_ -> "--exclude-intervals ${bed_}" }.join(" ") : ""
    def contig_ploidy = contig_ploidy_table ? "--contig-ploidy-priors ${contig_ploidy_table}" : ""
    def model = ploidy_model ? "--model ${ploidy_model}" : ""
    def input_list = counts.collect { count -> "--input ${count}" }.join(" ")

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK DetermineGermlineContigPloidy] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"
    export PYTENSOR_FLAGS="base_compiledir=\$PWD"
    export MPLCONFIGDIR="\$PWD"
    export XDG_CACHE_HOME="\$PWD"
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        DetermineGermlineContigPloidy \\
        ${input_list} \\
        --output ./ \\
        --output-prefix ${prefix} \\
        ${intervals} \\
        ${exclude} \\
        ${contig_ploidy} \\
        ${model} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-calls
    touch ${prefix}-model
    """
}

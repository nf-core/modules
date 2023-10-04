process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'
    
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "nf-core/universc:1.2.5.1"

    input:
    tuple val(meta), path(reads)
    path  reference


    output:
    tuple val(meta), path("sample-${meta.id}/outs/*"), emit: outs
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def sample_arg = meta.samples.unique().join(",")
    if ( reads instanceof Path || reads.size() != 2 ) {
        error "UNIVERSC module only supports 2 paired end read files. Input: ${reads}"
    }
    """
    universc \\
        --id ${meta.id} \\
        --read1 ${reads[0]} --read2 ${reads[1]} \\
        --reference ${reference} \\
        --jobmode "local" \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        $args

    echo !! > sample-${meta.id}/outs/_invocation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger:  \$(echo \$(cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g' ))
        universc:  \$(echo \$(bash /universc/launch_universc.sh --version | grep version | grep universc  | sed 's/^.* //g' ))
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    mkdir -p "${meta.id}/"
    mkdir -p "${meta.id}/SC_RNA_COUNTER_CS"
    mkdir -p "${meta.id}/journal"
    touch ${meta.id}/_log.txt
    touch ${meta.id}/test.mri.tgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger:  \$(echo \$(cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g' ))
        universc:  \$(echo \$(bash /universc/launch_universc.sh --version | grep version | grep universc | sed 's/^.* //g' ))
    END_VERSIONS
    """
}

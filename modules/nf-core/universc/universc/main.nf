process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/nfcore/universc:1.2.5.1"

    input:
    tuple val(meta), path(reads)
    path  reference

    output:
    tuple val(meta), path("${meta.id}", type: "dir"), emit: outs
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def sample_arg     = meta.samples.unique().join(",")
    def reference_name = reference.name
    // Older versions of this module supported a single read but Universc does not
    if ( reads instanceof Path || reads.size() != 2 ) {
        error "UNIVERSC module only supports paired end reads"
    }
    /*
    def input_reads    = "--read1 ${reads[0]} --read2 ${reads[1]}"
    */
    def technology_cmd = meta.technology ? "--technology ${meta.technology}" : ""
    def chemistry_cmd  = meta.chemistry  ? "--chemistry ${meta.chemistry}"   : ""
    """
    universc \\
        $args \\
        --id ${meta.id} \\
        --read1 ${reads[0]} --read2 ${reads[1]} \\
        ${technology_cmd} \\
        ${chemistry_cmd} \\
        --reference ${reference_name} \\
        --jobmode local \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger:  \$(echo \$(cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g' ))
        universc:  \$(echo \$(bash /universc/launch_universc.sh --version | grep version | grep universc  | sed 's/^.* //g' ))
    END_VERSIONS
    """


    stub:
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

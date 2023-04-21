process BCL2FASTQ {
    tag {"$meta.lane" ? "$meta.id"+"."+"$meta.lane" : "$meta.id" }
    label 'process_high'

    container "nfcore/bcl2fastq:2.20.0.422"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BCL2FASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("**[!Undetermined]_S*_R?_00?.fastq.gz"), emit: fastq
    tuple val(meta), path("**[!Undetermined]_S*_I?_00?.fastq.gz"), optional:true, emit: fastq_idx
    tuple val(meta), path("**Undetermined_S0*_R?_00?.fastq.gz")  , optional:true, emit: undetermined
    tuple val(meta), path("**Undetermined_S0*_I?_00?.fastq.gz")  , optional:true, emit: undetermined_idx
    tuple val(meta), path("Reports")                             , emit: reports
    tuple val(meta), path("Stats")                               , emit: stats
    tuple val(meta), path("InterOp/*.bin")                       , emit: interop
    path("versions.yml")                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def untar_dir = run_dir.toString().endsWith(".tar.gz") ? "tar -xzvf ${run_dir}" : ''
    def input_dir = untar_dir ? run_dir.toString() - '.tar.gz' : run_dir
    """
    $untar_dir

    bcl2fastq \\
        $args \\
        --output-dir . \\
        --runfolder-dir ${input_dir} \\
        --sample-sheet ${samplesheet} \\
        --processing-threads ${task.cpus}

    cp -r ${input_dir}/InterOp .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcl2fastq: \$(bcl2fastq -V 2>&1 | grep -m 1 bcl2fastq | sed 's/^.*bcl2fastq v//')
    END_VERSIONS
    """
}

process BCLCONVERT {
    tag {"$meta.lane" ? "$meta.id"+"."+"$meta.lane" : "$meta.id" }
    label 'process_high'

    container "nfcore/bclconvert:4.0.3"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BCLCONVERT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("**[!Undetermined]_S*_L00?_R?_00?.fastq.gz")  ,emit: fastq
    tuple val(meta), path("**_S*_L00?_I?_00?.fastq.gz")                 ,optional:true ,emit: fastq_idx
    tuple val(meta), path("Undetermined_S0_L00?_R?_00?.fastq.gz")       ,optional:true ,emit: undetermined
    tuple val(meta), path("Undetermined_S0_L00?_I?_00?.fastq.gz")       ,optional:true, emit: undetermined_idx
    tuple val(meta), path("Reports")                                    ,emit: reports
    tuple val(meta), path("Logs")                                       ,emit: logs
    tuple val(meta), path("**/InterOp/*.bin")                           ,emit: interop
    path("versions.yml")                                                ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcl-convert \\
        $args \\
        --output-directory . \\
        --bcl-input-directory ${run_dir} \\
        --sample-sheet ${samplesheet} \\
        --bcl-num-parallel-tiles ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(bcl-convert -V 2>&1 | head -n 1 | sed 's/^.*Version //')
    END_VERSIONS
    """
}

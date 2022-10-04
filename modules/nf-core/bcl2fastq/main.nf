process BCL2FASTQ {
    tag {"$meta.lane" ? "$meta.id"+"."+"$meta.lane" : "$meta.id" }
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using bcl2fastq. Please use docker or singularity containers."
    }
    container "nfcore/bcl2fastq:2.20.0.422"

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("**[!Undetermined]_S*_L00?_R?_00?.fastq.gz")  ,emit: fastq
    tuple val(meta), path("**_S*_L00?_I?_00?.fastq.gz")                 ,optional:true ,emit: fastq_idx
    tuple val(meta), path("Undetermined_S0_L00?_R?_00?.fastq.gz")       ,optional:true ,emit: undetermined
    tuple val(meta), path("Undetermined_S0_L00?_I?_00?.fastq.gz")       ,optional:true, emit: undetermined_idx
    tuple val(meta), path("Reports")                                    ,emit: reports
    tuple val(meta), path("Stats")                                      ,emit: stats
    tuple val(meta), path("**/InterOp/*.bin")                           ,emit: interop
    path("versions.yml")                                                ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcl2fastq \\
        $args \\
        --output-dir . \\
        --runfolder-dir ${run_dir} \\
        --sample-sheet ${samplesheet} \\
        --processing-threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcl2fastq: \$(bcl2fastq -V 2>&1 | grep -m 1 bcl2fastq | sed 's/^.*bcl2fastq v//')
    END_VERSIONS
    """
}

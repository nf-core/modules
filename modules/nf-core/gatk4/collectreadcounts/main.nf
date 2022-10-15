process GATK4_COLLECTREADCOUNTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    def container_image = "gatk4:4.2.6.1--hdfd78af_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("*.hdf5"), optional: true, emit: hdf5
    tuple val(meta), path("*.tsv") , optional: true, emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def reference = fasta ? "--reference $fasta" : ""
    def extension = args.contains("--format HDF5") ? "hdf5" :
                    args.contains("--format TSV")  ? "tsv" :
                    "hdf5"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK COLLECTREADCOUNTS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" CollectReadCounts \\
        --input $input \\
        --intervals $intervals \\
        --output ${prefix}.$extension \\
        $reference \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

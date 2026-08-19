process CELLRANGERARC_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/nf-core/cellranger-arc:2.0.2"

    input:
    tuple val(meta), val(sample_type), val(sub_sample), path(reads, stageAs: "fastqs/*")
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${prefix}/outs/**"), emit: outs
    tuple val(meta), path("${prefix}_lib.csv"), emit: lib
    tuple val("${task.process}"), val('cellrangerarc'), eval("cellranger-arc --version 2>&1 | sed 's/cellranger-arc cellranger-arc-//'"), emit: versions_cellrangerarc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def reference_name = reference.name
    def sample_types = sample_type.join(",")
    def sample_names = sub_sample.join(",")
    def lib_csv = meta.id + "_lib.csv"

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    fastq_folder=\$(readlink -f fastqs)

    python3 <<CODE

    sample_types = "${sample_types}".split(",")
    sample_names = "${sample_names}".split(",")
    unique_samples_names = set(sample_names)

    lib_csv = open("${lib_csv}", "w")
    lib_csv.write("fastqs,sample,library_type")

    for i in range(0, len(sample_types)):
        if sample_names[i] in unique_samples_names:
            unique_samples_names.remove(
                sample_names[i]
            )  # this has to be done to account for different Lane files (e.g., L002)
            if sample_types[i] == "gex":
                lib_csv.write("\\n{},{},{}".format("\${fastq_folder}", sample_names[i], "Gene Expression"))
            else:
                lib_csv.write("\\n{},{},{}".format("\${fastq_folder}", sample_names[i], "Chromatin Accessibility"))

    lib_csv.close()

    print("Wrote lib.csv file to {}".format("${lib_csv}"))

    CODE

    cellranger-arc \\
        count \\
        --id='${prefix}' \\
        --libraries=$lib_csv \\
        --reference=$reference_name \\
        --localcores=$task.cpus \\
        --localmem=${task.memory.toGiga()} \\
        $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "${prefix}/outs/"
    touch ${prefix}/outs/fake_file.txt
    touch ${prefix}_lib.csv
    """
}

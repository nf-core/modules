process SMNCOPYNUMBERCALLER {
    tag "$meta.id"
    label 'process_low'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used with SMNCopyNumberCaller at the moment. Please use Docker or Singularity containers."
    }
    container "clinicalgenomics/smncopynumbercaller:v1.1.2"

    input:
    tuple val(meta), path(bam)

    output:
    // QUESTION: what are good emit names?
    tuple val(meta), path ("*.txt"), emit: manifest_file
    tuple val(meta), path("*.tsv"), emit: smncopynumber_tsv
    tuple val(meta), path("*.json"), emit: run_metrics_json
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cpus = task.cpus
    def prefix = task.ext.prefix ?: "${meta.id}"
    def manifest_in = file(task.workDir+'/manifest.txt')
    manifest_in.text = bam.join("\n")
    def out_dir = task.workDir
    def genome_version = task.ext.genome_version
        // in the config file: def genome_version = genome.contains('GRch37') ? '37' : '38'

    """
    smn_caller.py \\
        $args \\
        --manifest $manifest_in \\
        --genome $genome_version \\
        --prefix $prefix \\
        --outDir $out_dir \\
        --threads $cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smncopynumbercaller: \$(echo \$(SMNCopyNumberCaller commit 3e67e3b on Feb 8, 2020 2>&1))
    END_VERSIONS
    """
}

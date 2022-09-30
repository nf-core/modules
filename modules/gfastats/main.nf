process GFASTATS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gfastats=1.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfastats:1.3.4--hd03093a_0':
        'quay.io/biocontainers/gfastats:1.3.4--hd03093a_0' }"

    input:
    tuple val(meta), path(assembly)   // input.[fasta|fastq|gfa][.gz]
    val genome_size                   // estimated genome size for NG* statistics (optional).
    val target                        // target specific sequence by header, optionally with coordinates (optional).
    path agp                          // -a --agp-to-path <file> converts input agp to path and replaces existing paths.
    path include_bed                  // -i --include-bed <file> generates output on a subset list of headers or coordinates in 0-based bed format.
    path exclude_bed                  // -e --exclude-bed <file> opposite of --include-bed. They can be combined (no coordinates).
    path instructions                 // -k --swiss-army-knife <file> set of instructions provided as an ordered list.

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def agpfile = agp ? "--agp-to-path $agp" : ""
    def ibed    = include_bed ? "--include-bed $include_bed" : ""
    def ebed    = exclude_bed ? "--exclude-bed $exclude_bed" : ""
    def manipulations = instructions ? "--swiss-army-knife $instructions" : ""
    """
    gfastats \\
        $args \\
        --threads $task.cpus \\
        $agpfile \\
        $ibed \\
        $ebed \\
        $manipulations \\
        $assembly \\
        $genome_size \\
        $target \\
        > ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """
}

process KMA_KMA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fc6c961562aef21c24b4f2330d9cd7e9bbda162b0d584a5cd5428e0b725e0d6/data':
        'community.wave.seqera.io/library/kma:1.5.0--eb093e0381fb59ea' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val (interleaved)

    output:
    tuple val(meta), path("*.res")    , optional: true, emit: res     // Results overview
    tuple val(meta), path("*.fsa")    , optional: true, emit: fsa     // Consensus sequences (disabled via '-nc')
    tuple val(meta), path("*.aln")    , optional: true, emit: aln     // Consensus alignments (disabled via '-na')
    tuple val(meta), path("*.frag.gz"), optional: true, emit: frag    // Read mapping information (disabled via '-nf')
    tuple val(meta), path("*.mat.gz") , optional: true, emit: matrix  // Base counts (only if -matrix is enabled)
    tuple val(meta), path("*.vcf.gz") , optional: true, emit: vcf
    tuple val(meta), path("*.sam")    , optional: true, emit: sam
    tuple val(meta), path("*.spa")    , optional: true, emit: spa
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Handle different read formats
    def read_command = interleaved ? "-int ${reads}" : meta.single_end ? "-i ${reads}" : "-ipe ${reads}"
    def sam_output = args.contains("-sam") ? "> ${prefix}.sam" : ''
    """
    # Determine index base name by looking for required index files
    INDEX_BASE=""
    INDEX_FILES=\$(ls ${index}/*.seq.b 2>/dev/null || ls ${index}/*/*.seq.b 2>/dev/null || echo "")

    if [ -n "\$INDEX_FILES" ]; then
        # Extract the base name from the first matching index file
        INDEX_FILE=\$(echo "\$INDEX_FILES" | head -n 1)
        INDEX_BASE=\${INDEX_FILE%.seq.b}
        echo "Using index base: \$INDEX_BASE"
    else
        # If no *.seq.b files found, try to check if the index itself is the base name
        if [ -f "${index}.seq.b" ]; then
            INDEX_BASE="${index}"
            echo "Using index base: \$INDEX_BASE"
        else
            echo "Error: Could not find proper KMA index files" >&2
            exit 1
        fi
    fi

    # FIXME: https://github.com/nf-core/modules/pull/7251
    # Run kma and handle exit code 95 as success
    kma \\
        ${read_command} \\
        -o ${prefix} \\
        -t_db \$INDEX_BASE \\
        $args \\
        $sam_output || [ \$? -eq 95 ]

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v 2>&1) | sed 's/^KMA-//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_alignments = args.contains('-Sparse') ?
        "touch ${prefix}.spa" :
        "touch ${prefix}.res; touch ${prefix}.fsa; touch ${prefix}.aln; echo \"\" | gzip > ${prefix}.frag.gz"
    def create_mat = args.contains('-mat') ? "echo \"\" | gzip > ${prefix}.mat.gz" : ""
    def create_vcf = args.contains('-vcf') ? "echo \"\" | gzip > ${prefix}.vcf.gz" : ""
    def create_sam = args.contains('-sam') ? "touch ${prefix}.sam" : ""
    """
    ${create_alignments}
    ${create_mat}
    ${create_vcf}
    ${create_sam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v 2>&1) | sed 's/^KMA-//')
    END_VERSIONS
    """
}

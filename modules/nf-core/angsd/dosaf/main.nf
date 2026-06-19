process ANGSD_DOSAF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--h13024bc_4':
        'quay.io/biocontainers/angsd:0.940--h13024bc_4' }"

    input:
    tuple val(meta), path(bams), path(bam_indices)
    tuple val(meta2), path(reference_fasta), path(reference_fai)
    tuple val(meta3), path(ancestral_fasta), path(ancestral_fai) // Optional. Provides ancestral state for unfolded SFS.
    tuple val(meta4), path(error_file) // Optional. Required for SYK model (-GL 4) only.
    tuple val(meta5), path(inbreeding_coefficients) // Optional. Required for -doSAF 2 (inbreeding-aware mode).

    output:
    tuple val(meta), path("*.saf.idx"), path("*.saf.pos.gz"), path("*.saf.gz"), emit: saf
    tuple val("${task.process}"), val('angsd'), eval("angsd 2>&1 | sed '1!d;s/.*version: //;s/ .*//'"), emit: versions_angsd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Validate GL model
    def gl_model = args.contains("-GL 1") ? 1 :
                   args.contains("-GL 2") ? 2 :
                   args.contains("-GL 3") ? 3 :
                   args.contains("-GL 4") ? 4 : 0
    if (gl_model == 0) {
        error(
            "ANGSD_DOSAF: No valid GL model selected. Please specify one of: -GL 1, -GL 2, -GL 3, or -GL 4 in ext.args."
            )
    }

    // Validate doSAF mode
    def dosaf_mode = args.contains("-doSAF 1") ? 1 :
                     args.contains("-doSAF 2") ? 2 : 0
    if (dosaf_mode == 0) {
        error(
            "ANGSD_DOSAF: No valid doSAF mode selected. Please specify -doSAF 1 or -doSAF 2 in ext.args."
            )
    }

    // Validate doSAF 2 has inbreeding coefficients
    if (dosaf_mode == 2 && !inbreeding_coefficients) {
        error(
            "ANGSD_DOSAF: -doSAF 2 requires inbreeding coefficients (indF). Please supply an inbreeding_coefficients file."
            )
    }

    // Validate GL 4 has error file
    if (gl_model == 4 && !error_file) {
        error(
            "ANGSD_DOSAF: -GL 4 (SYK model) requires an error file. Please supply an error_file."
            )
    }

    // Compute -minInd dynamically as a fraction of the population's BAM list size
    // Default set to 0 as this mimics angsd default for minInd
    def frac        = (task.ext.args2 ?: '0') as Double
    def n_samples   = bams instanceof List ? bams.size() : 1
    def min_ind_arg = frac > 0 ? "-minInd ${Math.ceil(n_samples * frac).toInteger()}" : ''

    // Touch fai indices to ensure they are newer than their fasta (ANGSD requirement)
    def touch_ref = reference_fai  ? "sleep 1 && touch ${reference_fai}"  : ''
    def touch_anc = ancestral_fai  ? "sleep 1 && touch ${ancestral_fai}"  : ''

    // Reference / ancestral args.
    // If ancestral_fasta supplied: use as -anc and pass reference as -ref.
    // If not supplied: use reference as -anc and -ref (suitable for folded SFS via realSFS -fold 1).
    def ref = ancestral_fasta ? "-anc ${ancestral_fasta} -ref ${reference_fasta}" : "-anc ${reference_fasta} -ref ${reference_fasta}"

    // Optional args
    def indF_arg   = inbreeding_coefficients ? "-indF ${inbreeding_coefficients}" : ""
    def errors_arg = error_file              ? "-errors ${error_file}"            : ""

    // Shared preamble
    def preamble = """
    ${touch_ref}
    ${touch_anc}
    printf '%s\\n' ${bams} > bamlist.txt
    """

    if (gl_model == 3) {
        // SOAPsnp requires a two-step process:
        // 1. Calibration step with -minQ 0 to generate the calibration matrix in angsd_tmpdir/
        // 2. doSAF step using the calibration matrix
        // -GL 3 is hardcoded in the calibration step to avoid passing all other args
        """
        ${preamble}

        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            -minQ 0 \\
            -GL 3 \\
            -ref ${reference_fasta} \\
            -out ${prefix}

        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            ${ref} \\
            ${indF_arg} \\
            ${args} \\
            ${min_ind_arg} \\
            -out ${prefix}
        """

    } else if (gl_model == 4) {
        // SYK model requires an error file and -doCounts 1
        """
        ${preamble}

        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            ${ref} \\
            ${indF_arg} \\
            ${errors_arg} \\
            -doCounts 1 \\
            ${args} \\
            ${min_ind_arg} \\
            -out ${prefix}
        """

    } else {
        // GL 1 (SAMtools) or GL 2 (GATK) - single angsd call
        """
        ${preamble}

        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            ${ref} \\
            ${indF_arg} \\
            ${args} \\
            ${min_ind_arg} \\
            -out ${prefix}
        """
    }

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.saf.idx
    echo "" | gzip > ${prefix}.saf.pos.gz
    echo "" | gzip > ${prefix}.saf.gz
    """
}

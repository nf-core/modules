process ANGSD_DOSAF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--h13024bc_4':
        'quay.io/biocontainers/angsd:0.940--h13024bc_4' }"

    input:
    tuple val(meta), path(bams), path(bam_indices)
    tuple path(reference_fasta), path(reference_fai)
    tuple path(ancestral_fasta), path(ancestral_fai) // Optional. Provides ancestral state for unfolded SFS.
    path(error_file) // Optional. Required for SYK model (-GL 4) only.
    path(inbreeding_coefficients) // Optional. Required for -doSAF 2 (inbreeding-aware mode).
    val(gl_model)
    val(dosaf_mode)

    output:
    tuple val(meta), path("*.saf.idx"), path("*.saf.pos.gz"), path("*.saf.gz"), emit: saf
    tuple val("${task.process}"), val('angsd'), eval("angsd 2>&1 | sed '1!d;s/.*version: //;s/ .*//'"), emit: versions_angsd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Validate GL model
    if (!(gl_model in [1, 2, 3, 4])) {
        error(
            "ANGSD_DOSAF: Invalid GL model '${gl_model}'. Please pass one of 1, 2, 3, or 4 as the gl_model input."
            )
    }

    // Validate doSAF mode
    if (!(dosaf_mode in [1, 2])) {
        error(
            "ANGSD_DOSAF: Invalid doSAF mode '${dosaf_mode}'. Please pass 1 or 2 as the dosaf_mode input."
            )
    }

    // Validate doSAF 2 has inbreeding coefficients
    if (dosaf_mode == 2 && !inbreeding_coefficients) {
        error(
            "ANGSD_DOSAF: -doSAF 2 requires inbreeding coefficients (indF). Please supply an inbreeding_coefficients file."
            )
    }

    // Validate doSAF 2 called alongside -doMajoMinor and -doMaf
    if (dosaf_mode == 2 && !(args.contains("-doMajorMinor"))) {
        error(
            "ANGSD_DOSAF: -doSAF 2 requires -doMajorMinor to be specified. Please include in ext.args."
            )
    }
    if (dosaf_mode == 2 && !(args.contains("-doMaf"))) {
        error(
            "ANGSD_DOSAF: -doSAF 2 requires -doMaf to be specified. Please include in ext.args."
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
    def ref_anc_arg = ancestral_fasta ? "-anc ${ancestral_fasta} -ref ${reference_fasta}" : "-anc ${reference_fasta} -ref ${reference_fasta}"

    // Optional args
    def indF_arg   = inbreeding_coefficients ? "-indF ${inbreeding_coefficients}" : ""
    def errors_arg = error_file              ? "-errors ${error_file}"            : ""

    // Shared preamble
    def preamble = """
    ${touch_ref}
    ${touch_anc}
    printf '%s\\n' ${bams} > bamlist.txt
    """

// SOAPsnp (-GL 3) needs a calibration pass first to generate the matrix
// angsd reads back in during the doSAF step. -GL 3 is hardcoded here
// since no other args should reach the calibration call.
def calibration = ''
if (gl_model == 3) {
    calibration = """
    angsd \\
        -nThreads ${task.cpus} \\
        -bam bamlist.txt \\
        -minQ 0 \\
        -GL ${gl_model} \\
        -ref ${reference_fasta} \\
        -out ${prefix}
    """
}

// SYK model (-GL 4) needs an error file and -doCounts 1
def gl4_args = gl_model == 4 ? "${errors_arg} -doCounts 1" : ''

"""
${preamble}
${calibration}
angsd \\
    -nThreads ${task.cpus} \\
    -bam bamlist.txt \\
    -GL ${gl_model} \\
    -doSAF ${dosaf_mode} \\
    ${ref} \\
    ${indF_arg} \\
    ${gl4_args} \\
    ${args} \\
    ${min_ind_arg} \\
    -out ${prefix}
"""

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.saf.idx
    echo "" | gzip > ${prefix}.saf.pos.gz
    echo "" | gzip > ${prefix}.saf.gz
    """
}

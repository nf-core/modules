process ANGSD_GL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--hce60e53_2':
        'biocontainers/angsd:0.940--hce60e53_2' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(fasta)      //Optionally
    tuple val(meta3), path(error_file) //Optionally. Used for SYK model only.

    output:
    tuple val(meta), path("*.{glf,beagle}.gz"), emit: genotype_likelihood
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def GL_model = args.contains("-GL 1") ? 1 : args.contains("-GL 2") ? 2 : args.contains("-GL 3") ? 3 : args.contains("-GL 4") ? 4 : 0
    def ref = fasta ? "-ref ${fasta}" : ''                     // Use reference fasta if provided
    def errors = error_file ? "-errors ${error_file}" : ''     // Only applies to SYK model
    def output_mode = args.contains("-doGlf") ? "" : '-doGlf 1' // Default to outputting binary glf (10 log likelihoods) if not set in args
    // NOTE: GL is specified within args, so is not provided as a separate argument

    if (GL_model != 3 && GL_model != 4) {
        """
        ls -1 *.bam > bamlist.txt

        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            $args \\
            $ref \\
            $output_mode \\
            -out ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            angsd: \$(echo \$(angsd 2>&1) | grep 'angsd version' | head -n 1 | sed 's/.*version: //g;s/ .*//g')
        END_VERSIONS
        """
    } else if (GL_model == 3) {
        // No args for this part.
        // GL is hardcoded to 3 here to avoid passing all other arguments to the calibration step
        """
        ls -1 *.bam > bamlist.txt

        ## SOAPsnp model
        ## First get the calibration matrix. minQ MUST be 0 for this step. Will create the directory angsd_tmpdir/ with the required files for the next step.
        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            -minQ 0 \\
            -GL 3 \\
            $ref \\
            -out ${prefix}

        ## Then run the model
        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            $args \\
            $ref \\
            $output_mode \\
            -out ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            angsd: \$(echo \$(angsd 2>&1) | grep 'angsd version' | head -n 1 | sed 's/.*version: //g;s/ .*//g')
        END_VERSIONS
        """
    } else if (GL_model == 4) {
        """
        ls -1 *.bam > bamlist.txt

        ## SYK model
        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            $args \\
            $ref \\
            $output_mode \\
            $errors \\
            -doCounts 1 \\
            -out ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            angsd: \$(echo \$(angsd 2>&1) | grep 'angsd version' | head -n 1 | sed 's/.*version: //g;s/ .*//g')
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.glf
    gzip ${prefix}.glf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep 'angsd version' | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}

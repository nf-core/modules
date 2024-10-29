def is_compressed(filename, extensions = ['gz', 'gzip'] as Set<String>) {
    return extensions.contains(filename.extension)
}

def is_sequence(filename, extensions = ['fasta', 'fa', 'fastq', 'fq'] as Set<String>) {
    return extensions.contains(filename.extension)
}

def get_simple_name(filename) {
    // If the input has a gzip compression, strip the file extension.
    if (is_compressed(filename)) {
        filename = java.nio.file.Paths.get(filename.baseName)
    }
    // Strip an expected sequencing file extension.
    if (is_sequence(filename)) {
        filename = java.nio.file.Paths.get(filename.baseName)
    } else {
        throw new Exception("Unrecognized sequencing file extension '${filename.extension}'.")
    }
    return filename.name
}

def get_single_prefix(sequence_file) {
    return get_simple_name(sequence_file)
}

def get_paired_prefix(first, _second) {
    // Paired sequencing files usually have a suffix of `_1` or `.1` in their name,
    // which we remove.
    return get_simple_name(first) - ~/[._]1$/
}


process KRAKENUNIQ_PRELOADEDKRAKENUNIQ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.4--pl5321h6dccd9a_2':
        'biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_2' }"

    input:
    tuple val(meta), path(sequences)
    val sequence_type
    path db
    val ram_chunk_size
    val save_output_reads
    val report_file
    val save_output

    output:
    tuple val(meta), path("*.classified.${sequence_type}.gz")  , optional:true, emit: classified_reads
    tuple val(meta), path("*.unclassified.${sequence_type}.gz"), optional:true, emit: unclassified_reads
    tuple val(meta), path('*.krakenuniq.classified.txt')       , optional:true, emit: classified_assignment
    tuple val(meta), path('*.krakenuniq.report.txt')           , emit: report
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assert sequence_type in ['fasta', 'fastq']

    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''

    classified   = meta.single_end ? "\${PREFIX}.classified.${sequence_type}"   : "\${PREFIX}.merged.classified.${sequence_type}"
    unclassified = meta.single_end ? "\${PREFIX}.unclassified.${sequence_type}" : "\${PREFIX}.merged.unclassified.${sequence_type}"
    classified_option = save_output_reads ? "--classified-out \"${classified}\"" : ''
    unclassified_option = save_output_reads ? "--unclassified-out \"${unclassified}\"" : ''
    def output_option = save_output ? '--output "\${PREFIX}.krakenuniq.classified.txt"' : ''
    def report = report_file ? '--report-file "\${PREFIX}.krakenuniq.report.txt"' : ''
    compress_reads_command = save_output_reads ? "find . -name '*.${sequence_type}' -print0 | xargs -0 -t -P ${task.cpus} -I % gzip --no-name %" : ''

    if (meta.single_end) {
        command_input = sequences.collect { seq -> "${seq} ${get_single_prefix(seq)}" }
        """
        cat <<-END_OPTIONS > options.txt
        ${command_input.join('\n')}
        END_OPTIONS

        krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        while IFS= read -r SEQ PREFIX; do
            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                $args2 \\
                "\${SEQ}"
        done

        $compress_reads_command

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        command_input = sequences.collect { first, second -> "${first} ${second} ${get_paired_prefix(first, second)}" }
        """
        cat <<-END_OPTIONS > options.txt
        ${command_input.join('\n')}
        END_OPTIONS

        krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        while IFS= read -r SEQ_1 SEQ_2 PREFIX; do
            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                --paired \\
                $args2 \\
                "\${SEQ_1}" "\${SEQ_2}"
        done

        $compress_reads_command

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    }

    stub:
    assert sequence_type in ['fasta', 'fastq']

    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''

    classified   = meta.single_end ? "\${PREFIX}.classified.${sequence_type}"   : "\${PREFIX}.merged.classified.${sequence_type}"
    unclassified = meta.single_end ? "\${PREFIX}.unclassified.${sequence_type}" : "\${PREFIX}.merged.unclassified.${sequence_type}"
    classified_option = save_output_reads ? "--classified-out \"${classified}\"" : ''
    unclassified_option = save_output_reads ? "--unclassified-out \"${unclassified}\"" : ''
    def output_option = save_output ? '--output "\${PREFIX}.krakenuniq.classified.txt"' : ''
    def report = report_file ? '--report-file "\${PREFIX}.krakenuniq.report.txt"' : ''
    compress_reads_command = save_output_reads ? "find . -name '*.${sequence_type}' -print0 | xargs -0 -t -P ${task.cpus} -I % gzip --no-name %" : ''

    if (meta.single_end) {
        command_input = sequences.collect { seq -> "${seq} ${get_single_prefix(seq)}" }
        """
        cat <<-END_OPTIONS > options.txt
        ${command_input.join('\n')}
        END_OPTIONS

        echo krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        create_file() {
            echo '<3 nf-core' > "\$1"
        }

        create_gzip_file() {
            echo '<3 nf-core' | gzip -n > "\$1"
        }

        while IFS= read -r SEQ PREFIX; do
            echo "\${SEQ}"
            echo "\${PREFIX}"

            echo krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                $args2 \\
                "\${SEQ}"

            create_file "\${PREFIX}.krakenuniq.classified.txt"
            create_file "\${PREFIX}.krakenuniq.report.txt"
            create_gzip_file "\${PREFIX}.classified.${sequence_type}.gz"
            create_gzip_file "\${PREFIX}.unclassified.${sequence_type}.gz"
        done

        echo "$compress_reads_command"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        command_input = sequences.collect { first, second -> "${first} ${second} ${get_paired_prefix(first, second)}" }
        """
        cat <<-END_OPTIONS > options.txt
        ${command_input.join('\n')}
        END_OPTIONS

        echo krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        create_file() {
            echo '<3 nf-core' > "\$1"
        }

        create_gzip_file() {
            echo '<3 nf-core' | gzip -n > "\$1"
        }

        while IFS= read -r SEQ_1 SEQ_2 PREFIX; do
            echo "\${SEQ_1} \${SEQ_2}"
            echo "\${PREFIX}"

            echo krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                \${OPTIONS} \\
                --paired \\
                $args2 \\
                "\${SEQ_1}" "\${SEQ_2}"

            create_file "\${PREFIX}.krakenuniq.classified.txt"
            create_file "\${PREFIX}.krakenuniq.report.txt"
            create_gzip_file "\${PREFIX}.merged.classified.${sequence_type}.gz"
            create_gzip_file "\${PREFIX}.merged.unclassified.${sequence_type}.gz"
        done

        echo "$compress_reads_command"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    }
}

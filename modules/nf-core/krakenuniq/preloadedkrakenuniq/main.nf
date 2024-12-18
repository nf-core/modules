process KRAKENUNIQ_PRELOADEDKRAKENUNIQ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.4--pl5321h6dccd9a_2':
        'biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_2' }"

    input:
    // We stage sequencing files in a sub-directory so we don't accidentally gzip them later.
    tuple val(meta), path(sequences, name: 'sequences/*'), val(prefixes)
    val sequence_type
    path db
    val ram_chunk_size
    val save_output_reads
    val report_file
    val save_output

    output:
    tuple val(meta), path("*.classified.${sequence_type}.gz")  , optional:true, emit: classified_reads
    tuple val(meta), path("*.unclassified.${sequence_type}.gz"), optional:true, emit: unclassified_reads
    tuple val(meta), path("*.krakenuniq.classified.txt")       , optional:true, emit: classified_assignment
    tuple val(meta), path("*.krakenuniq.report.txt")           , emit: report
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assert sequence_type in ['fasta', 'fastq']
    sequences = sequences instanceof List ? sequences : [sequences]

    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''

    classified   = meta.single_end ? "\${PREFIX}.classified.${sequence_type}"   : "\${PREFIX}.merged.classified.${sequence_type}"
    unclassified = meta.single_end ? "\${PREFIX}.unclassified.${sequence_type}" : "\${PREFIX}.merged.unclassified.${sequence_type}"
    classified_option = save_output_reads ? "--classified-out \"${classified}\"" : ''
    unclassified_option = save_output_reads ? "--unclassified-out \"${unclassified}\"" : ''
    def output_option = save_output ? '--output "\${PREFIX}.krakenuniq.classified.txt"' : ''
    def report = report_file ? '--report-file "\${PREFIX}.krakenuniq.report.txt"' : ''
    compress_reads_command = save_output_reads ? "find . -maxdepth 0 -name '*.${sequence_type}' -print0 | xargs -0 -t -P ${task.cpus} -I % gzip --no-name %" : ''
    def command_inputs_file = '.inputs.txt'

    if (meta.single_end) {
        assert sequences.size() == prefixes.size()
        command_inputs = [sequences, prefixes].transpose().collect { seq, prefix -> "${seq}\t${prefix}" }

        """
        # Store the batch of samples for later command input.
        cat <<-END_INPUTS > ${command_inputs_file}
        ${command_inputs.join('\n        ')}
        END_INPUTS

        # Preload the KrakenUniq database into memory.
        krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        # Run the KrakenUniq classification on each sample in the batch.
        while IFS='\t' read -r SEQ PREFIX; do
            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                $args2 \\
                "\${SEQ}"
        done < ${command_inputs_file}

        $compress_reads_command

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        assert sequences.size() / 2 == prefixes.size()
        command_inputs = [sequences.collate(2), prefixes].transpose().collect { pair, prefix -> "${pair[0]}\t${pair[1]}\t${prefix}" }

        """
        # Store the batch of samples for later command input.
        cat <<-END_INPUTS > ${command_inputs_file}
        ${command_inputs.join('\n        ')}
        END_INPUTS

        # Preload the KrakenUniq database into memory.
        krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        # Run the KrakenUniq classification on each sample in the batch.
        while IFS='\t' read -r FIRST_SEQ SECOND_SEQ PREFIX; do
            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                --paired \\
                $args2 \\
                "\${FIRST_SEQ}" "\${SECOND_SEQ}"
        done < ${command_inputs_file}

        $compress_reads_command

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    }

    stub:
    assert sequence_type in ['fasta', 'fastq']
    sequences = sequences instanceof List ? sequences : [sequences]

    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''

    classified   = meta.single_end ? "\${PREFIX}.classified.${sequence_type}"   : "\${PREFIX}.merged.classified.${sequence_type}"
    unclassified = meta.single_end ? "\${PREFIX}.unclassified.${sequence_type}" : "\${PREFIX}.merged.unclassified.${sequence_type}"
    classified_option = save_output_reads ? "--classified-out \"${classified}\"" : ''
    unclassified_option = save_output_reads ? "--unclassified-out \"${unclassified}\"" : ''
    def output_option = save_output ? '--output "\${PREFIX}.krakenuniq.classified.txt"' : ''
    def report = report_file ? '--report-file "\${PREFIX}.krakenuniq.report.txt"' : ''
    compress_reads_command = save_output_reads ? "find . -name '*.${sequence_type}' -print0 | xargs -0 -t -P ${task.cpus} -I % gzip --no-name %" : ''
    def command_inputs_file = '.inputs.txt'

    if (meta.single_end) {
        assert sequences.size() == prefixes.size()
        command_inputs = [sequences, prefixes].transpose().collect { seq, prefix -> "${seq}\t${prefix}" }

        """
        # Store the batch of samples for later command input.
        cat <<-END_INPUTS > ${command_inputs_file}
        ${command_inputs.join('\n        ')}
        END_INPUTS

        # Preload the KrakenUniq database into memory.
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

        # Run the KrakenUniq classification on each sample in the batch.
        while IFS='\t' read -r SEQ PREFIX; do
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
        done < ${command_inputs_file}

        echo "$compress_reads_command"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        assert sequences.size() / 2 == prefixes.size()
        command_inputs = [sequences.collate(2), prefixes].transpose().collect { pair, prefix -> "${pair[0]}\t${pair[1]}\t${prefix}" }

        """
        # Store the batch of samples for later command input.
        cat <<-END_INPUTS > ${command_inputs_file}
        ${command_inputs.join('\n        ')}
        END_INPUTS

        # Preload the KrakenUniq database into memory.
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

        # Run the KrakenUniq classification on each sample in the batch.
        while IFS='\t' read -r FIRST_SEQ SECOND_SEQ PREFIX; do
            echo krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                --paired \\
                $args2 \\
                "\${FIRST_SEQ}" "\${SECOND_SEQ}"

            create_file "\${PREFIX}.krakenuniq.classified.txt"
            create_file "\${PREFIX}.krakenuniq.report.txt"
            create_gzip_file "\${PREFIX}.merged.classified.${sequence_type}.gz"
            create_gzip_file "\${PREFIX}.merged.unclassified.${sequence_type}.gz"
        done < ${command_inputs_file}

        echo "$compress_reads_command"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    }
}

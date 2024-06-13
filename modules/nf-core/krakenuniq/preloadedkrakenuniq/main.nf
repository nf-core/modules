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
        """
        krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        strip_suffix() {
            local result=\$1
            # Strip any file extensions.
            echo "\${result%%.*}"
        }

        printf "%s\\n" ${sequences} | while read FASTQ; do \\
            PREFIX="\$(strip_suffix "\${FASTQ}")"

            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                $args2 \\
                "\${FASTQ}"
        done

        $compress_reads_command

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        strip_suffix() {
            local result
            read result
            # Strip any trailing dot or underscore.
            result="\${result%_}"
            echo "\${result%.}"
        }

        printf "%s %s\\n" ${sequences} | while read FASTQ; do \\
            read -r -a FASTQ <<< "\${FASTQ}"
            PREFIX="\$(printf "%s\\n" "\${FASTQ[@]}" |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' | strip_suffix)"

            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                --paired \\
                $args2 \\
                "\${FASTQ[@]}"
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
        """
        echo krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        strip_suffix() {
            local result=\$1
            # Strip any file extensions.
            echo "\${result%%.*}"
        }

        create_file() {
            echo '<3 nf-core' > "\$1"
        }

        create_gzip_file() {
            echo '<3 nf-core' | gzip -n > "\$1"
        }

        printf "%s\\n" ${sequences} | while read FASTQ; do \\
            echo "\${FASTQ}"
            PREFIX="\$(strip_suffix "\${FASTQ}")"
            echo "\${PREFIX}"

            echo krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                $args2 \\
                "\${FASTQ}"

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
        """
        echo krakenuniq \\
            $args \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus

        strip_suffix() {
            local result
            read result
            # Strip any trailing dot or underscore.
            result="\${result%_}"
            echo "\${result%.}"
        }

        create_file() {
            echo '<3 nf-core' > "\$1"
        }

        create_gzip_file() {
            echo '<3 nf-core' | gzip -n > "\$1"
        }

        printf "%s %s\\n" ${sequences} | while read FASTQ; do \\
            read -r -a FASTQ <<< "\${FASTQ}"
            echo "\${FASTQ[@]}"
            PREFIX="\$(printf "%s\\n" "\${FASTQ[@]}" |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' | strip_suffix)"
            echo "\${PREFIX}"

            echo krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                --paired \\
                $args2 \\
                "\${FASTQ[@]}"

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

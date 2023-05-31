process KRAKENUNIQ_PRELOADEDKRAKENUNIQ {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::krakenuniq=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.2--pl5321h19e8d03_0':
        'biocontainers/krakenuniq:1.0.2--pl5321h19e8d03_0' }"

    input:
    tuple val(meta), path(fastqs)
    path  db
    val ram_chunk_size
    val save_output_fastqs
    val report_file
    val save_output

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classified.txt')        , optional:true, emit: classified_assignment
    tuple val(meta), path('*report.txt')                           , emit: report

    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''

    def classified   = meta.single_end ? '"\${PREFIX}.classified.fastq"'   : '"\${PREFIX}.classified#.fastq"'
    def unclassified = meta.single_end ? '"\${PREFIX}.unclassified.fastq"' : '"\${PREFIX}.unclassified#.fastq"'
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ''
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ''
    def output_option = save_output ? '--output "\${PREFIX}.krakenuniq.classified.txt"' : ''
    def report = report_file ? '--report-file "\${PREFIX}.krakenuniq.report.txt"' : ''
    def compress_reads_command = save_output_fastqs ? 'gzip --no-name *.fastq' : ''
    if (meta.single_end) {
        """
        krakenuniq \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus \\
            $args

        strip_suffix() {
            local result=\$1
            # Strip any file extensions.
            echo "\${result%%.*}"
        }

        printf "%s\\n" ${fastqs} | while read FASTQ; do \\
            PREFIX="\$(strip_suffix "\${FASTQ}")"

            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                $output_option \\
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
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus \\
            $args

        strip_suffix() {
            local result
            read result
            # Strip any trailing dot or underscore.
            result="\${result%_}"
            echo "\${result%.}"
        }

        printf "%s %s\\n" ${fastqs} | while read FASTQ; do \\
            read -r -a FASTQ <<< "\${FASTQ}"
            PREFIX="\$(printf "%s\\n" "\${FASTQ[@]}" |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' | strip_suffix)"

            krakenuniq \\
                --db $db \\
                --threads $task.cpus \\
                $report \\
                $output_option \\
                $unclassified_option \\
                $classified_option \\
                $output_option \\
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
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''

    def classified   = meta.single_end ? '"\${PREFIX}.classified.fastq"'   : '"\${PREFIX}.classified#.fastq"'
    def unclassified = meta.single_end ? '"\${PREFIX}.unclassified.fastq"' : '"\${PREFIX}.unclassified#.fastq"'
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ''
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ''
    def output_option = save_output ? '--output "\${PREFIX}.krakenuniq.classified.txt"' : ''
    def report = report_file ? '--report-file "\${PREFIX}.krakenuniq.report.txt"' : ''
    def compress_reads_command = save_output_fastqs ? 'gzip --no-name *.fastq' : ''
    if (meta.single_end) {
        """
        echo krakenuniq \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus \\
            $args

        strip_suffix() {
            local result=\$1
            # Strip any file extensions.
            echo "\${result%%.*}"
        }

        printf "%s\\n" ${fastqs} | while read FASTQ; do \\
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
                $output_option \\
                $args2 \\
                "\${FASTQ}"

            touch "\${PREFIX}.classified.fastq.gz"
            touch "\${PREFIX}.krakenuniq.classified.txt"
            touch "\${PREFIX}.krakenuniq.report.txt"
            touch "\${PREFIX}.unclassified.fastq.gz"
        done

        echo $compress_reads_command

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        echo krakenuniq \\
            --db $db \\
            --preload \\
            --preload-size $ram_chunk_size \\
            --threads $task.cpus \\
            $args

        strip_suffix() {
            local result
            read result
            # Strip any trailing dot or underscore.
            result="\${result%_}"
            echo "\${result%.}"
        }

        printf "%s %s\\n" ${fastqs} | while read FASTQ; do \\
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
                $output_option \\
                --paired \\
                $args2 \\
                "\${FASTQ[@]}"

            touch "\${PREFIX}.classified_1.fastq.gz" "\${PREFIX}.classified_2.fastq.gz"
            touch "\${PREFIX}.krakenuniq.classified.txt"
            touch "\${PREFIX}.krakenuniq.report.txt"
            touch "\${PREFIX}.unclassified_1.fastq.gz" "\${PREFIX}.unclassified_2.fastq.gz"
        done

        echo $compress_reads_command

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
        END_VERSIONS
        """
    }
}

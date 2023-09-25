process NANOCOMP {
    label 'process_medium'

    conda "bioconda:nanocomp=1.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp:1.21.0--pyhdfd78af_0':
        'biocontainers/nanocomp:1.21.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(filelist)

    output:
    tuple val(meta), path("*NanoComp-report.html"), emit: report_html
    tuple val(meta), path("*NanoComp_lengths_violin.html"), emit: lengths_violin_html
    tuple val(meta), path("*NanoComp_log_length_violin.html"), emit: log_length_violin_html
    tuple val(meta), path("*NanoComp_N50.html"), emit: n50_html
    tuple val(meta), path("*NanoComp_number_of_reads.html"), emit: number_of_reads_html
    tuple val(meta), path("*NanoComp_OverlayHistogram.html"), emit: overlay_histogram_html
    tuple val(meta), path("*NanoComp_OverlayHistogram_Normalized.html"), emit: overlay_histogram_normalized_html
    tuple val(meta), path("*NanoComp_OverlayLogHistogram.html"), emit: overlay_log_histogram_html
    tuple val(meta), path("*NanoComp_OverlayLogHistogram_Normalized.html"), emit: overlay_log_histogram_normalized_html
    tuple val(meta), path("*NanoComp_total_throughput.html"), emit: total_throughput_html
    tuple val(meta), path("*NanoComp_quals_violin.html"), emit: quals_violin_html, optional: true
    tuple val(meta), path("*NanoComp_OverlayHistogram_Identity.html"), emit: overlay_histogram_identity_html, optional: true
    tuple val(meta), path("*NanoComp_OverlayHistogram_PhredScore.html"), emit: overlay_histogram_phredscore_html, optional: true
    tuple val(meta), path("*NanoComp_percentIdentity_violin.html"), emit: percent_identity_violin_html, optional: true
    tuple val(meta), path("*NanoComp_ActivePoresOverTime.html"), emit: active_pores_over_time_html, optional: true
    tuple val(meta), path("*NanoComp_CumulativeYieldPlot_Gigabases.html"), emit: cumulative_yield_plot_gigabases_html, optional: true
    tuple val(meta), path("*NanoComp_sequencing_speed_over_time.html"), emit: sequencing_speed_over_time_html, optional: true
    tuple val(meta), path("*NanoStats.txt"), emit: stats_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (prefix == ""){
        prefixflag = ""
    } else {
        prefixflag = "--prefix " + prefix
    }

    //determine input file type
    filetypes = []
    for (file in filelist){
        tokenized_filename = file.getName().tokenize('.')
        if (tokenized_filename.size() < 2){
            throw new java.lang.IndexOutOfBoundsException("Every input file to nanocomp has to have a file ending.")
        }

        first_namepart = true
        extension_found = false

        for (namepart in tokenized_filename){
            if (namepart == ""){
                continue
            }

            // prevent the file name to be seen as extension
            if (first_namepart == true){
                first_namepart = false
                continue
            }

            if (["fq","fastq"].contains(namepart)){
                filetypes.add("fastq")
                extension_found = true
                break
            } else if (["fasta", "fna", "ffn", "faa", "frn", "fa"].contains(namepart)) {
                filetypes.add("fasta")
                extension_found = true
                break
            } else if (namepart == "bam") {
                filetypes.add("bam")
                extension_found = true
                break
            } else if (namepart == "txt") {
                filetypes.add("summary")
                extension_found = true
                break
            }
        }

        if (extension_found == false){
            throw new java.lang.IllegalArgumentException("There was no suitable filetype found for " + file.getName() +
            ". NanoComp only accepts fasta (fasta, fna, ffn, faa, frn, fa), fastq (fastq, fq), bam and Nanopore sequencing summary (txt).")
        }
    }

    filetypes.unique()
    if (filetypes.size() < 1){
        throw new java.lang.IllegalArgumentException("There was no suitable filetype found in NanoComp input. Please use fasta, fastq, bam or Nanopore sequencing summary.")
    }
    if (filetypes.size() > 1){
        throw new java.lang.IllegalArgumentException("You gave different filetypes to NanoComp. Please use only *one* of fasta, fastq, bam or Nanopore sequencing summary.")
    }
    filetype = filetypes[0]

    """
    NanoComp \\
        --$filetype $filelist \\
        --threads $task.cpus \\
        $prefixflag \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocomp: \$(echo \$(NanoComp --version 2>&1) | sed 's/^.*NanoComp //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch versions.yml
    touch "${prefix}"NanoComp_lengths_violin.html
    touch "${prefix}"NanoComp_log_length_violin.html
    touch "${prefix}"NanoComp_N50.html
    touch "${prefix}"NanoComp_number_of_reads.html
    touch "${prefix}"NanoComp_OverlayHistogram.html
    touch "${prefix}"NanoComp_OverlayHistogram_Normalized.html
    touch "${prefix}"NanoComp_OverlayLogHistogram.html
    touch "${prefix}"NanoComp_OverlayLogHistogram_Normalized.html
    touch "${prefix}"NanoComp-report.html
    touch "${prefix}"NanoComp_total_throughput.html
    touch "${prefix}"NanoComp_quals_violin.html
    touch "${prefix}"NanoComp_OverlayHistogram_Identity.html
    touch "${prefix}"NanoComp_OverlayHistogram_PhredScore.html
    touch "${prefix}"NanoComp_percentIdentity_violin.html
    touch "${prefix}"NanoComp_ActivePoresOverTime.html
    touch "${prefix}"NanoComp_CumulativeYieldPlot_Gigabases.html
    touch "${prefix}"NanoComp_sequencing_speed_over_time.html
    touch "${prefix}"NanoStats.txt
    """
}

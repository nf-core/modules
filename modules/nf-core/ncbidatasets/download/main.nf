process NCBIDATASETS_DOWNLOAD {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "conda-forge::ncbi-datasets-cli=15.11.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-datasets-cli:15.11.0':
        'staphb/ncbi-datasets:15.11.0' }"

    input:

    tuple val(meta), val(extra_args)

    output:
    // Required output
    tuple val(meta), path("*.zip")   ,             emit: zip
    path "versions.yml"              ,             emit: versions
    // Optional outputs, depends on tool usage
    tuple val(meta), path("cds.fna.gz"),           emit: cds,        optional: true
    tuple val(meta), path("gene.fna.gz"),          emit: gene,       optional: true
    tuple val(meta), path("*genomic.fna.gz"),      emit: genome,     optional: true
    tuple val(meta), path("genomic.gbff.gz"),      emit: gbff,       optional: true
    tuple val(meta), path("genomic.gff.gz"),       emit: gff,        optional: true
    tuple val(meta), path("genomic.gtf.gz"),       emit: gtf,        optional: true
    tuple val(meta), path("protein.faa.gz"),       emit: protein,    optional: true
    tuple val(meta), path("rna.fna.gz"),           emit: rna,        optional: true
    tuple val(meta), path("3p_utr.fna.gz"),        emit: utr_3p,     optional: true
    tuple val(meta), path("5p_utr.fna.gz"),        emit: utr_5p,     optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def valid_commands = ["genome accession", "genome taxon", "gene accession", "gene gene-id", "gene symbol", "gene taxon", "virus genome taxon", "virus genome accession", "virus protein"]
    if (!valid_commands.contains(meta.command)) {
        error "Unsupported command: ${meta.command} "
    }


    def args = task.ext.args ?: ''
    def args_from_csv = extra_args ?: ''
    args_from_csv = args_from_csv.replaceAll("^\"|\"\$", "")
    def prefix = task.ext.prefix ?: "${meta.id.replaceAll(' ', '_')}"

    // In case of downloading based on taxon identifiers:
    // multiple genomes/proteomes etc may exist in different files in the ncbi .zip file
    // Therefore, only emit .zip file and do not emit extracted files
    if (meta.command =~ /taxon/)
        """

        [ -e /usr/local/ssl/cacert.pem ] && export SSL_CERT_FILE=/usr/local/ssl/cacert.pem

        datasets download \\
            --filename ${prefix}_ncbi_dataset.zip \\
            ${args} \\
            ${args_from_csv} \\
            ${meta.command} "${meta.id}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ncbi-datasets-cli: \$(echo \$(datasets --version 2>&1) | sed 's/datasets version: //' )
        END_VERSIONS
        """

    // Download .zip file and emit
    // Also, extract (optional) faa, fna, gff, gtf and gbff files and emit them in their respective channels
    else
        """
        [ -e /usr/local/ssl/cacert.pem ] && export SSL_CERT_FILE=/usr/local/ssl/cacert.pem

        datasets download \\
            --filename ${prefix}_ncbi_dataset.zip \\
            ${args} \\
            ${args_from_csv} \\
            ${meta.command} ${meta.id}


        # extract all .fna and .faa files to current folder and compress them
        unzip -j ${prefix}_ncbi_dataset.zip
        find \\( -name "*.faa" -o -name "*.fna" -o -name "*.gff" -o -name "*.gtf"  -o -name "*.gbff" \\) -exec gzip -n {} \\;

        # rename cds_from_genomic
        if [ -f cds_from_genomic.fna.gz ]; then
            mv cds_from_genomic.fna.gz cds.fna.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ncbi-datasets-cli: \$(echo \$(datasets --version 2>&1) | sed 's/datasets version: //' )
        END_VERSIONS
        """



    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

process MASURCA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "ecoflowucl/masurca:v4.1.4"

    input:
    tuple val(meta), path(illumina), path(jump), path(pacbio), path(nanopore)
    val fragment_mean
    val fragment_stdev
    val jump_mean
    val jump_stdev
    val extend_jump_reads
    val graph_kmer_size
    val use_linking_mates
    val lhe_coverage
    val mega_reads_one_pass
    val limit_jump_coverage
    val ca_parameters
    val close_gaps
    val jf_size

    output:
    tuple val(meta), path("assemble.sh"), emit: script
    tuple val(meta), path("*scaffolds.fa.gz"), emit: scaffolds
    tuple val(meta), path("*_masurca_config.txt"), emit: config
    tuple val(meta), path("*-masurca.log"), emit: log
    tuple val("${task.process}"), val('masurca'), eval("masurca --version | sed 's/version //g'"), topic: versions, emit: versions_masurca

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    //get input reads with absolute paths - illumina are mandatory, jump/pacbio/nanopore are optional
    def illumina_reads = [illumina].flatten().collect { it.toRealPath() }.join(' ')
    def jump_reads = jump ? [jump].flatten().collect { it.toRealPath() }.join(' ') : ""
    def pacbio_file = pacbio ? pacbio.toRealPath() : ""
    def nanopore_file = nanopore ? nanopore.toRealPath() : ""

    // Build the config file
    def config_lines = []

    // DATA section
    config_lines << "DATA"
    config_lines << "#Illumina paired end reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>"
    config_lines << "#if single-end, do not specify <reverse_reads>"
    config_lines << "#MUST HAVE Illumina paired end reads to use MaSuRCA"
    config_lines << "PE= pe ${fragment_mean} ${fragment_stdev} ${illumina_reads}"

    // Jump/mate pair reads (optional)
    if (jump_reads) {
        config_lines << "#Illumina mate pair reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>"
        config_lines << "JUMP= sh ${jump_mean} ${jump_stdev} ${jump_reads}"
    }


    // PacBio and Nanopore reads handling
    def long_reads_concat = ""
    if (pacbio_file && nanopore_file) {
        config_lines << "#if you have both PacBio and Nanopore, supply both as NANOPORE type"
        long_reads_concat = "${prefix}_long_reads.fastq.gz"
        config_lines << "NANOPORE= ${long_reads_concat}"
    } else if (pacbio_file) {
        config_lines << "#PacBio/CCS reads must be in a single fasta or fastq file with absolute path"
        config_lines << "PACBIO=${pacbio_file}"
    } else if (nanopore_file) {
        config_lines << "#Nanopore reads must be in a single fasta or fastq file with absolute path"
        config_lines << "NANOPORE=${nanopore_file}"
    }

    config_lines << "END"
    config_lines << ""

    // PARAMETERS section
    config_lines << "PARAMETERS"
    config_lines << "#set this to 1 if your Illumina jumping library reads are shorter than 100bp"
    config_lines << "EXTEND_JUMP_READS=${extend_jump_reads}"
    config_lines << "#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content"
    config_lines << "GRAPH_KMER_SIZE = ${graph_kmer_size}"
    config_lines << "#set this to 1 for all Illumina-only assemblies"
    config_lines << "#set this to 0 if you have more than 15x coverage by long reads (Pacbio or Nanopore) or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc)"
    config_lines << "USE_LINKING_MATES = ${use_linking_mates}"
    config_lines << "#use at most this much coverage by the longest Pacbio or Nanopore reads, discard the rest of the reads"
    config_lines << "#can increase this to 30 or 35 if your reads are short (N50<7000bp)"
    config_lines << "LHE_COVERAGE=${lhe_coverage}"
    config_lines << "#set to 0 (default) to do two passes of mega-reads for slower, but higher quality assembly, otherwise set to 1"
    config_lines << "MEGA_READS_ONE_PASS=${mega_reads_one_pass}"
    config_lines << "#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms"
    config_lines << "LIMIT_JUMP_COVERAGE = ${limit_jump_coverage}"
    config_lines << "#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically."
    config_lines << "#CABOG ASSEMBLY ONLY: set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms."
    config_lines << "CA_PARAMETERS = ${ca_parameters}"
    config_lines << "#CABOG ASSEMBLY ONLY: whether to attempt to close gaps in scaffolds with Illumina  or long read data"
    config_lines << "CLOSE_GAPS=${close_gaps}"
    config_lines << "#number of cpus to use, set this to the number of CPUs/threads per node you will be using"
    config_lines << "NUM_THREADS = ${task.cpus}"
    config_lines << "#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20"
    config_lines << "JF_SIZE = ${jf_size}"
    config_lines << "#ILLUMINA ONLY. Set this to 1 to use SOAPdenovo contigging/scaffolding module."
    config_lines << "#Assembly will be worse but will run faster. Useful for very large (>=8Gbp) genomes from Illumina-only data"
    config_lines << "SOAP_ASSEMBLY=0"
    config_lines << "#If you are doing Hybrid Illumina paired end + Nanopore/PacBio assembly ONLY (no Illumina mate pairs or OTHER frg files)."
    config_lines << "#Set this to 1 to use Flye assembler for final assembly of corrected mega-reads."
    config_lines << "#A lot faster than CABOG, AND QUALITY IS THE SAME OR BETTER."
    config_lines << "#Works well even when MEGA_READS_ONE_PASS is set to 1."
    config_lines << "#DO NOT use if you have less than 15x coverage by long reads."
    config_lines << "FLYE_ASSEMBLY=0"
    config_lines << "END"

    def config_content = config_lines.join('\n')

    """
    # Write the config file
cat > ${prefix}_masurca_config.txt <<-CONFIG_EOF
${config_content}
CONFIG_EOF

    # Concatenate long reads if both PacBio and Nanopore are present
    ${long_reads_concat ? "cat ${pacbio_file} ${nanopore_file} > ${long_reads_concat}" : ""}

    # Generate assembly script
    masurca ${prefix}_masurca_config.txt

    ./assemble.sh > ${prefix}-masurca.log 2>&1

    if [ -f CA*/primary.genome.scf.fasta ]; then
        gzip -cn CA*/primary.genome.scf.fasta > ${prefix}.scaffolds.fa.gz
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p CA
    touch assemble.sh
    touch ${prefix}_masurca_config.txt
    echo | gzip > ${prefix}.scaffolds.fa.gz
    touch ${prefix}-masurca.log
    """
}

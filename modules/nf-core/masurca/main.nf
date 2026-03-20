// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process MASURCA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cf/cf6402ed20c3b089ab88cd8884ddace90693501453a515f9188ae681e8ca8556/data':
        'community.wave.seqera.io/library/masurca:4.1.4--d05ef74c4881d55c' }"

    input:
    tuple val(meta), path(illumina), path(jump), path(pacbio), path(nanopore), path(other_reads), path(reference_genome)
    val fragment_mean
    val fragment_stdev
    val jump_mean
    val jump_stdev


    output:
    tuple val(meta), path("${prefix}/assemble.sh")                               , emit: script
    tuple val(meta), path("${prefix}/CA*/final.genome.scf.fasta"), optional: true, emit: scaffolds
    tuple val(meta), path("${prefix}/CA*/final.genome.ctg.fasta"), optional: true, emit: contigs
    tuple val(meta), path("${prefix}/flye/assembly.fasta")       , optional: true, emit: flye_assembly
    tuple val(meta), path("${prefix}/*_masurca_config.txt")                      , emit: config
    tuple val("${task.process}"), val('masurca'), eval("masurca --version | sed 's/version //g'"), topic: versions, emit: versions_masurca

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    //get input reads with absolute paths - illumina are mandatory, jump/pacbio/nanopore are optional
    def illumina_reads = meta.single_end ? "\$(readlink -f ${illumina})" : "\$(readlink -f ${illumina[0]}) \$(readlink -f ${illumina[1]})"
    def jump_reads = jump ? "\$(readlink -f ${jump[0]}) \$(readlink -f ${jump[1]})" : ""
    def pacbio_file = pacbio ? "\$(readlink -f ${pacbio})" : ""
    def nanopore_file = nanopore ? "\$(readlink -f ${nanopore})" : ""
    def other_reads_file = other_reads ? "\$(readlink -f ${other_reads})" : ""
    def reference_genome_file = reference_genome ? "\$(readlink -f ${reference_genome})" : ""

    // Configuration parameters with defaults from task.ext
    def extend_jump_reads = task.ext.extend_jump_reads != null ? task.ext.extend_jump_reads : 0
    def graph_kmer_size = task.ext.graph_kmer_size ?: 'auto'
    def use_linking_mates = task.ext.use_linking_mates != null ? task.ext.use_linking_mates : 0
    def lhe_coverage = task.ext.lhe_coverage ?: 25
    def mega_reads_one_pass = task.ext.mega_reads_one_pass != null ? task.ext.mega_reads_one_pass : 0
    def limit_jump_coverage = task.ext.limit_jump_coverage ?: 300
    def ca_parameters = task.ext.ca_parameters ?: 'cgwErrorRate=0.15'
    def close_gaps = task.ext.close_gaps != null ? task.ext.close_gaps : 1
    def jf_size = task.ext.jf_size ?: 200000000
    def soap_assembly = task.ext.soap_assembly != null ? task.ext.soap_assembly : 0
    def flye_assembly = task.ext.flye_assembly != null ? task.ext.flye_assembly : 0
    """
    echo "DATA" > ${prefix}_masurca_config.txt
    echo "#Illumina paired end reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>" >> ${prefix}_masurca_config.txt
    echo "#if single-end, do not specify <reverse_reads>" >> ${prefix}_masurca_config.txt
    echo "#MUST HAVE Illumina paired end reads to use MaSuRCA" >> ${prefix}_masurca_config.txt
    echo "PE= pe ${fragment_mean} ${fragment_stdev} ${illumina_reads}" >> ${prefix}_masurca_config.txt

    # Jump/mate pair reads (optional)
    if [ -n "${jump_reads}" ]; then
        echo "#Illumina mate pair reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>" >> ${prefix}_masurca_config.txt
        echo "JUMP= sh ${jump_mean} ${jump_stdev} ${jump_reads}" >> ${prefix}_masurca_config.txt
    fi

    # PacBio and Nanopore reads handling
    # If both exist, concatenate them and supply as NANOPORE (per MaSuRCA docs)
    if [ -n "${pacbio_file}" ] && [ -n "${nanopore_file}" ]; then
        echo "#if you have both PacBio and Nanopore, supply both as NANOPORE type" >> ${prefix}_masurca_config.txt
        cat ${pacbio_file} ${nanopore_file} > ${prefix}_long_reads.fastq.gz
        echo "NANOPORE=\$(readlink -f ${prefix}_long_reads.fastq.gz)" >> ${prefix}_masurca_config.txt
    elif [ -n "${pacbio_file}" ]; then
        echo "#PacBio/CCS reads must be in a single fasta or fastq file with absolute path" >> ${prefix}_masurca_config.txt
        echo "PACBIO=${pacbio_file}" >> ${prefix}_masurca_config.txt
    elif [ -n "${nanopore_file}" ]; then
        echo "#Nanopore reads must be in a single fasta or fastq file with absolute path" >> ${prefix}_masurca_config.txt
        echo "NANOPORE=${nanopore_file}" >> ${prefix}_masurca_config.txt
    fi

    # Other reads (optional) - Sanger, 454, etc.
    if [ -n "${other_reads_file}" ]; then
        echo "#Other reads (Sanger, 454, etc) one frg file, concatenate your frg files into one if you have many" >> ${prefix}_masurca_config.txt
        echo "OTHER=${other_reads_file}" >> ${prefix}_masurca_config.txt
    fi

    # Reference genome (optional) - for synteny-assisted assembly
    if [ -n "${reference_genome_file}" ]; then
        echo "#synteny-assisted assembly, concatenate all reference genomes into one reference.fa; works for Illumina-only data" >> ${prefix}_masurca_config.txt
        echo "REFERENCE=${reference_genome_file}" >> ${prefix}_masurca_config.txt
    fi

    echo "END" >> ${prefix}_masurca_config.txt


    echo "" >> ${prefix}_masurca_config.txt
    echo "PARAMETERS" >> ${prefix}_masurca_config.txt
    echo "#set this to 1 if your Illumina jumping library reads are shorter than 100bp" >> ${prefix}_masurca_config.txt
    echo "EXTEND_JUMP_READS=${extend_jump_reads}" >> ${prefix}_masurca_config.txt
    echo "#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content" >> ${prefix}_masurca_config.txt
    echo "GRAPH_KMER_SIZE = ${graph_kmer_size}" >> ${prefix}_masurca_config.txt
    echo "#set this to 1 for all Illumina-only assemblies" >> ${prefix}_masurca_config.txt
    echo "#set this to 0 if you have more than 15x coverage by long reads (Pacbio or Nanopore) or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc)" >> ${prefix}_masurca_config.txt
    echo "USE_LINKING_MATES = ${use_linking_mates}" >> ${prefix}_masurca_config.txt
    echo "#use at most this much coverage by the longest Pacbio or Nanopore reads, discard the rest of the reads" >> ${prefix}_masurca_config.txt
    echo "#can increase this to 30 or 35 if your reads are short (N50<7000bp)" >> ${prefix}_masurca_config.txt
    echo "LHE_COVERAGE=${lhe_coverage}" >> ${prefix}_masurca_config.txt
    echo "#set to 0 (default) to do two passes of mega-reads for slower, but higher quality assembly, otherwise set to 1" >> ${prefix}_masurca_config.txt
    echo "MEGA_READS_ONE_PASS=${mega_reads_one_pass}" >> ${prefix}_masurca_config.txt
    echo "#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms" >> ${prefix}_masurca_config.txt
    echo "LIMIT_JUMP_COVERAGE = ${limit_jump_coverage}" >> ${prefix}_masurca_config.txt
    echo "#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically." >> ${prefix}_masurca_config.txt
    echo "#CABOG ASSEMBLY ONLY: set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms." >> ${prefix}_masurca_config.txt
    echo "CA_PARAMETERS = ${ca_parameters}" >> ${prefix}_masurca_config.txt
    echo "#CABOG ASSEMBLY ONLY: whether to attempt to close gaps in scaffolds with Illumina  or long read data" >> ${prefix}_masurca_config.txt
    echo "CLOSE_GAPS=${close_gaps}" >> ${prefix}_masurca_config.txt
    echo "#number of cpus to use, set this to the number of CPUs/threads per node you will be using" >> ${prefix}_masurca_config.txt
    echo "NUM_THREADS = ${task.cpus}" >> ${prefix}_masurca_config.txt
    echo "#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20" >> ${prefix}_masurca_config.txt
    echo "JF_SIZE = ${jf_size}" >> ${prefix}_masurca_config.txt
    echo "#ILLUMINA ONLY. Set this to 1 to use SOAPdenovo contigging/scaffolding module." >> ${prefix}_masurca_config.txt
    echo "#Assembly will be worse but will run faster. Useful for very large (>=8Gbp) genomes from Illumina-only data" >> ${prefix}_masurca_config.txt
    echo "SOAP_ASSEMBLY=${soap_assembly}" >> ${prefix}_masurca_config.txt
    echo "#If you are doing Hybrid Illumina paired end + Nanopore/PacBio assembly ONLY (no Illumina mate pairs or OTHER frg files)." >> ${prefix}_masurca_config.txt
    echo "#Set this to 1 to use Flye assembler for final assembly of corrected mega-reads." >> ${prefix}_masurca_config.txt
    echo "#A lot faster than CABOG, AND QUALITY IS THE SAME OR BETTER." >> ${prefix}_masurca_config.txt
    echo "#Works well even when MEGA_READS_ONE_PASS is set to 1." >> ${prefix}_masurca_config.txt
    echo "#DO NOT use if you have less than 15x coverage by long reads." >> ${prefix}_masurca_config.txt
    echo "FLYE_ASSEMBLY=${flye_assembly}" >> ${prefix}_masurca_config.txt
    echo "END" >> ${prefix}_masurca_config.txt
    
    # Generate assembly script
    masurca ${prefix}_masurca_config.txt
    
    # Create output directory and move files
    mkdir -p ${prefix}
    mv assemble.sh ${prefix}/
    mv ${prefix}_masurca_config.txt ${prefix}/
    chmod +x ${prefix}/assemble.sh
    
    # Run the assembly
    cd ${prefix}
    ./assemble.sh
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    mkdir -p ${prefix}/CA
    mkdir -p ${prefix}/flye
    touch ${prefix}/assemble.sh
    touch ${prefix}/${prefix}_masurca_config.txt
    touch ${prefix}/CA/final.genome.scf.fasta
    touch ${prefix}/CA/final.genome.ctg.fasta
    touch ${prefix}/flye/assembly.fasta
    """
}

//
// Runs FGBIO tools to remove UMI tags from FASTQ reads
// Convert them to unmapped BAM file, map them to the reference genome,
// use the mapped information to group UMIs and generate consensus reads
//

params.fastqtobam_options      = [:]
params.bam2fq_options          = [args: "-T RX"]
params.bwamem_options          = [args: "-p -C -M"]
params.samblaster_options      = [args: "-M --addMateTags", suffix:'_processed']
params.groupreadsbyumi_options = [:]
params.groupumistrategy        = "Adjacency"
params.umiconsensus_options    = [args: '-M 1 -S Coordinate', suffix: '_umiconsensus']

include { FGBIO_FASTQTOBAM                  as FASTQTOBAM }         from '../../../modules/fgbio/fastqtobam/main'                        addParams( options: params.fastqtobam_options )
include { SAMTOOLS_BAM2FQ                   as BAM2FASTQ }          from '../../../modules/samtools/bam2fq/main.nf'                   addParams( options: params.bam2fq_options )
include { BWA_INDEX }                                               from '../../../modules/bwa/index/main.nf'                         addParams( options: [:] )
include { BWA_MEM }                                                 from '../../../modules/bwa/mem/main'                                 addParams( options: params.bwamem_options )
include { SAMBLASTER }                                              from '../../../modules/samblaster/main'                              addParams( options: params.samblaster_options )
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI }    from '../../../modules/fgbio/groupreadsbyumi/main'                   addParams( options: params.groupreadsbyumi_options )
include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLUMICONSENSUS }   from '../../../modules/fgbio/callmolecularconsensusreads/main.nf' addParams( options: params.umiconsensus_options )


worfklow CREATE_UMI_CONSENSUS {
    take:
    reads                     // channel: [ val(meta), [ reads ] ]
    fasta                     // channel: /path/to/reference/fasta
    read_structure            // channel: val(read_structure)

    main:
    ch_versions = Channel.empty()
    fastqtobam_input = Channel.from(reads)


    // using information in val(read_structure) FASTQ reads are converted into
    // a tagged unmapped BAM file (uBAM)
    FASTQTOBAM ( fastqtobam_input, read_structure )
    ch_versions = ch_versions.mix(FASTQTOBAM.out.versions)

    // reference is indexed
    BWA_INDEX ( fasta )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    // in order to map uBAM using BWA MEM, we need to convert uBAM to FASTQ
    // but keep the appropriate UMI tags in the FASTQ comment field and produce
    // an interleaved FASQT file (hence, split = false)
    split = false
    BAM2FASTQ ( FASTQTOBAM.out.umibam, split )
    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)

    // appropriately tagged interleaved FASTQ reads are mapped to the reference
    BWA_MEM ( BAM2FASTQ.out.reads, BWA_INDEX.out.index )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    // samblaster is used in order to tag mates information in the BAM file
    // this is used in order to group reads by UMI
    SAMBLASTER ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMBLASTER.out.versions)

    // appropriately tagged reads are now grouped by UMI information
    GROUPREADSBYUMI ( SAMBLASTER.out.bam, params.groupumistrategy )
    ch_versions = ch_versions.mix(GROUPREADSBYUMI.out.versions)

    // using the above created groups, a consensus across reads in the same grou
    // can be called
    // this will emit a consensus BAM file
    CALLUMICONSENSUS ( GROUPREADSBYUMI.out.bam )
    ch_versions = ch_versions.mix(CALLUMICONSENSUS.out.versions)

    emit:
    ubam           = FASTQTOBAM.out.umibam          // channel: [ val(meta), [ bam ] ]
    groupbam       = GROUPREADSBYUMI.out.bam        // channel: [ val(meta), [ bam ] ]
    consensusbam   = CALLUMICONSENSUS.out.bam       // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}


//
// Runs FGBIO tools to remove UMI tags from FASTQ reads
// Convert them to unmapped BAM file, map them to the reference genome,
// use the mapped information to group UMIs and generate consensus reads
//

include { BWA_INDEX                         as BWAMEM1_INDEX       } from '../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM                           as BWAMEM1_MEM_PRE     } from '../../../modules/nf-core/bwa/mem/main.nf'
include { BWA_MEM                           as BWAMEM1_MEM_POST    } from '../../../modules/nf-core/bwa/mem/main.nf'
include { BWAMEM2_INDEX                                            } from '../../../modules/nf-core/bwamem2/index/main.nf'
include { BWAMEM2_MEM                       as BWAMEM2_MEM_PRE     } from '../../../modules/nf-core/bwamem2/mem/main.nf'
include { BWAMEM2_MEM                       as BWAMEM2_MEM_POST    } from '../../../modules/nf-core/bwamem2/mem/main.nf'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLUMICONSENSUS    } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_CALLDUPLEXCONSENSUSREADS    as CALLDUPLEXCONSENSUS } from '../../../modules/nf-core/fgbio/callduplexconsensusreads/main.nf'
include { FGBIO_FASTQTOBAM                  as FASTQTOBAM          } from '../../../modules/nf-core/fgbio/fastqtobam/main.nf'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUS     } from '../../../modules/nf-core/fgbio/filterconsensusreads/main.nf'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI     } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main.nf'
include { FGBIO_ZIPPERBAMS                  as ZIPPERBAMS_PRE      } from '../../../modules/nf-core/fgbio/zipperbams/main.nf'
include { FGBIO_ZIPPERBAMS                  as ZIPPERBAMS_POST     } from '../../../modules/nf-core/fgbio/zipperbams/main.nf'
include { SAMTOOLS_FASTQ                    as BAM2FASTQ_PRE       } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { SAMTOOLS_FASTQ                    as BAM2FASTQ_POST      } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { SAMTOOLS_SORT                     as SORTBAM             } from '../../../modules/nf-core/samtools/sort/main.nf'
workflow FASTQ_CREATE_UMI_CONSENSUS_FGBIO {

    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] [ val(meta), /path/to/reference/fasta ]
    bwa_index                 // channel: [optional]  [ val(meta), /path/to/reference/bwaindex ]
    dict                      // channel: [mandatory] [ val(meta), /path/to/reference/dictionary ]
    groupreadsbyumi_strategy  // string:  [mandatory] grouping strategy - default: "Adjacency"
    aligner                   // string:  [mandatory] "bwa-mem" or "bwa-mem2"
    duplex                    // bool:    [mandatory] true or false depending on UMI structure
    min_reads                 //          [mandatory] One integer (for non-duplex) or a string of up-to three space-separated numbers for duplex sequencing
    min_baseq                 // integer: [mandatory]
    max_base_error_rate       // integer: [mandatory] Maximum base error rate for consensus building

    main:

    ch_versions = Channel.empty()

    // reference is indexed if index not available in iGenomes - this is set in modules configuration
    // NB: this should exist in main workflow in a form like:
    // params.bwaindex = WorkflowMain.getGenomeAttribute(params, 'bwa')

    // using information in val(read_structure) FASTQ reads are converted into
    // a tagged unmapped BAM file (uBAM)
    // if the UMIs are present in read names instead of inline sequences
    // please make sure you adjust your config to include --extract-umis-from-read-names with ext.args
    // of the following step
    FASTQTOBAM ( reads )
    ch_versions = ch_versions.mix(FASTQTOBAM.out.versions)

    // in order to map uBAM using BWA MEM, we need to convert uBAM to FASTQ
    BAM2FASTQ_PRE ( FASTQTOBAM.out.bam, false )
    ch_versions = ch_versions.mix(BAM2FASTQ_PRE.out.versions)

    // the user can choose here to use either bwa-mem (default) or bwa-mem2
    aligned_bam = Channel.empty()

    if (aligner == "bwa-mem") {

        if(!bwa_index){
            BWAMEM1_INDEX ( fasta )
            ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
        }

        // sets bwaindex to correct input
        bwaindex    = bwa_index ?: BWAMEM1_INDEX.out.index
        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        // the aligner should be set with the following parameters "-p -K 150000000 -Y"
        // to be configured in ext.args of your config
        BWAMEM1_MEM_PRE ( BAM2FASTQ_PRE.out.fastq, bwaindex, fasta, false )
        ch_versions = ch_versions.mix(BWAMEM1_MEM_PRE.out.versions)
        aligned_bam = aligned_bam.mix(BWAMEM1_MEM_PRE.out.bam)
    } else {

        if(!bwa_index){
            BWAMEM2_INDEX ( fasta )
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }

        // sets bwaindex to correct input
        bwaindex    = bwa_index ?: BWAMEM2_INDEX.out.index
        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        // the aligner should be set with the following parameters "-p -K 150000000 -Y"
        // to be configured in ext.args of your config
        BWAMEM2_MEM_PRE ( BAM2FASTQ_PRE.out.fastq, bwaindex, fasta, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM_PRE.out.versions)
        aligned_bam = BWAMEM2_MEM_PRE.out.bam
    }

    // in order to tag mates information in the BAM file
    // FGBIO tool ZipperBams is used to merge info from mapped and unmapped BAM files
    ZIPPERBAMS_PRE ( FASTQTOBAM.out.bam, aligned_bam, fasta, dict )
    ch_versions = ch_versions.mix(ZIPPERBAMS_PRE.out.versions)

    // appropriately tagged reads are now grouped by UMI information
    // note that in tests ext.args has been set to recommended --edits 1
    // if UMIs are significantly longer (e.g. 20bp) or have more errors, --edits can be increased
    // explicit parameter strategy is Adjacency
    // For multiplex PCR and similar data where reads' genomic positions are fixed by the primers
    // it is recommended to use --strategy Identity to reduce runtime at the expense of lower accuracy
    // For duplex UMIs reads MUST be grouped using --strategy paired
    GROUPREADSBYUMI ( ZIPPERBAMS_PRE.out.bam, groupreadsbyumi_strategy )
    ch_versions = ch_versions.mix(GROUPREADSBYUMI.out.versions)

    // prepare output channel independently on UMI structure
    consensus_bam = Channel.empty()

    if (duplex){
        // this is executed if the library contains duplex UMIs
        CALLDUPLEXCONSENSUS ( GROUPREADSBYUMI.out.bam, min_reads, min_baseq )
        ch_versions = ch_versions.mix(CALLDUPLEXCONSENSUS.out.versions)
        consensus_bam =  CALLDUPLEXCONSENSUS.out.bam

    } else {
        // using the above created groups, a consensus across reads in the same group
        // can be called
        // this will emit a consensus BAM file
        CALLUMICONSENSUS ( GROUPREADSBYUMI.out.bam, min_reads, min_baseq )
        ch_versions = ch_versions.mix(CALLUMICONSENSUS.out.versions)
        consensus_bam =  CALLUMICONSENSUS.out.bam
    }

    FILTERCONSENSUS ( consensus_bam, fasta, min_reads, min_baseq, max_base_error_rate )
    ch_versions = ch_versions.mix(FILTERCONSENSUS.out.versions)

    // now the consensus uBAM needs to be converted into FASTQ again
    // to be aligned
    BAM2FASTQ_POST ( FILTERCONSENSUS.out.bam, false )
    ch_versions = ch_versions.mix(BAM2FASTQ_POST.out.versions)

    if (aligner == "bwa-mem") {
        // index made available through previous steps
        BWAMEM1_MEM_POST ( BAM2FASTQ_POST.out.fastq, bwaindex, fasta, false )
        ch_versions = ch_versions.mix(BWAMEM1_MEM_POST.out.versions)
        aligned_bam_post = BWAMEM1_MEM_POST.out.bam
    } else {
        // index made available through previous steps
        BWAMEM2_MEM_POST ( BAM2FASTQ_POST.out.fastq, bwaindex, fasta, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM_POST.out.versions)
        aligned_bam_post = BWAMEM2_MEM_POST.out.bam
    }

    // in order to tag mates information in the BAM file
    // FGBIO tool ZipperBams is used to merge info from mapped and unmapped BAM files
    ZIPPERBAMS_POST ( consensus_bam, aligned_bam_post, fasta, dict )
    ch_versions = ch_versions.mix(ZIPPERBAMS_POST.out.versions)

    // finally sort bam file
    SORTBAM ( ZIPPERBAMS_POST.out.bam, fasta )
    ch_versions = ch_versions.mix(SORTBAM.out.versions)

    emit:
    ubam               = FASTQTOBAM.out.bam             // channel: [ val(meta), [ bam ] ]
    groupbam           = GROUPREADSBYUMI.out.bam        // channel: [ val(meta), [ bam ] ]
    consensusbam       = consensus_bam                  // channel: [ val(meta), [ bam ] ]
    mappedconsensusbam = SORTBAM.out.bam                // channel: [ val(meta), [ bam ] ]
    versions           = ch_versions                    // channel: [ versions.yml ]
}

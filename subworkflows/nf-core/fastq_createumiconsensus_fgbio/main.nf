//
// Runs FGBIO tools to remove UMI tags from FASTQ reads
// Convert them to unmapped BAM file, map them to the reference genome,
// use the mapped information to group UMIs and generate consensus reads
//

include { BWAMEM2_INDEX                                            } from '../../../modules/nf-core/bwamem2/index/main.nf'
include { BWAMEM2_MEM                                              } from '../../../modules/nf-core/bwamem2/mem/main.nf'
include { BWA_INDEX                         as BWAMEM1_INDEX       } from '../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM                           as BWAMEM1_MEM         } from '../../../modules/nf-core/bwa/mem/main.nf'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLUMICONSENSUS    } from '../../../modules/nf-core/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_CALLDUPLEXCONSENSUSREADS    as CALLDUPLEXCONSENSUS } from '../../../modules/nf-core/fgbio/callduplexconsensusreads/main.nf'
include { FGBIO_FASTQTOBAM                  as FASTQTOBAM          } from '../../../modules/nf-core/fgbio/fastqtobam/main.nf'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUS     } from '../../../modules/nf-core/fgbio/filterconsensusreads/main.nf'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI     } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main.nf'
include { FGBIO_ZIPPERBAMS                  as ZIPPERBAMS          } from '../../../modules/nf-core/fgbio/zipperbams/main.nf'
include { SAMTOOLS_FASTQ                    as BAM2FASTQ           } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { SAMTOOLS_SORT                     as SORTBAM             } from '../../../modules/nf-core/samtools/sort/main.nf'


workflow FASTQ_CREATEUMICONSENSUS_FGBIO {

    take:
    reads                     // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                     // channel: [mandatory] /path/to/reference/fasta
    dict                      // channel: [mandatory] /path/to/reference/dictionary
    read_structure            // string:  [mandatory] "read_structure"
    groupreadsbyumi_strategy  // string:  [mandatory] grouping strategy - default: "Adjacency"
    aligner                   // string:  [mandatory] "bwa-mem" or "bwa-mem2"
    duplex                    // bool:    [mandatory] true or false depending on UMI structure

    main:

    ch_versions = Channel.empty()

    // using information in val(read_structure) FASTQ reads are converted into
    // a tagged unmapped BAM file (uBAM)
    FASTQTOBAM ( reads, read_structure )
    ch_versions = ch_versions.mix(FASTQTOBAM.out.version)

    // in order to map uBAM using BWA MEM, we need to convert uBAM to FASTQ
    BAM2FASTQ ( FASTQTOBAM.out.umibam )
    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)

    // the user can choose here to use either bwa-mem (default) or bwa-mem2
    aligned_bam = Channel.empty()

    if (aligner == "bwa-mem") {
        // reference is indexed if index not available in iGenomes - this is set in modules configuration
        // NB: this should exist in main workflow in a form like:
        // params.bwaindex = WorkflowMain.getGenomeAttribute(params, 'bwa')
        BWAMEM1_INDEX ( fasta )
        ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)

        // sets bwaindex to correct input
        bwaindex    = params.fasta ? params.bwaindex ? Channel.fromPath(params.bwaindex).collect() : BWAMEM1_INDEX.out.index : []
        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        BWAMEM1_MEM ( BAM2FASTQ.out.reads, bwaindex, false )
        ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions)
        aligned_bam = BWAMEM1_MEM.out.bam
    } else {
        // reference is indexed if index not available in iGenomes - this is set in modules configuration
        // NB: this should exist in main workflow in a form like:
        // params.bwaindex = WorkflowMain.getGenomeAttribute(params, 'bwa')
        BWAMEM2_INDEX ( fasta )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

        // sets bwaindex to correct input
        bwaindex    = params.fasta ? params.bwaindex ? Channel.fromPath(params.bwaindex).collect() : BWAMEM1_INDEX.out.index : []

        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        BWAMEM2_MEM ( BAM2FASTQ.out.reads, bwaindex, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        aligned_bam = BWAMEM2_MEM.out.bam
    }

    // in order to tag mates information in the BAM file
    // FGBIO tool ZipperBams is used to merge info from mapped and unmapped BAM files
    ZIPPERBAMS ( FASTQTOBAM.out.umibam, aligned_bam, fasta, dict )
    ch_versions = ch_versions.mix(ZIPPERBAMS.out.versions)

    // appropriately tagged reads are now grouped by UMI information
    GROUPREADSBYUMI ( ZIPPERBAMS.out.bam, groupreadsbyumi_strategy )
    ch_versions = ch_versions.mix(GROUPREADSBYUMI.out.versions)

    // please note in FILTERCONSENSUSREADS:
    // --min-reads is a required argument with no default
    // --min-base-quality is a required argument with no default
    // make sure they are specified via ext.args in your config
    FILTERCONSENSUS ( GROUPREADSBYUMI.out.bam, fasta )
    ch_versions = ch_versions.mix(FILTERCONSENSUS.out.versions)

    // prepare output channel independently on UMI structure
    consensus_bam = Channel.empty()

    if (duplex){
        // this is executed if the library contains duplex UMIs
        CALLDUPLEXCONSENSUS ( )


    } else {
        // using the above created groups, a consensus across reads in the same group
        // can be called
        // this will emit a consensus BAM file
        CALLUMICONSENSUS ( FILTERCONSENSUS.out.bam )
        ch_versions = ch_versions.mix(CALLUMICONSENSUS.out.versions)
        consensus_bam =  CALLUMICONSENSUS.out.bam
    }





    emit:
    ubam           = FASTQTOBAM.out.umibam          // channel: [ val(meta), [ bam ] ]
    groupbam       = GROUPREADSBYUMI.out.bam        // channel: [ val(meta), [ bam ] ]
    consensusbam   = CALLUMICONSENSUS.out.bam       // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}


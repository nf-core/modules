//
// Runs FGBIO tools to remove UMI tags from FASTQ reads
// Convert them to unmapped BAM file, map them to the reference genome,
// use the mapped information to group UMIs and generate consensus reads
//

include { BWAMEM2_INDEX                     } from '../../../modules/bwamem2/index/main.nf'
include { BWAMEM2_MEM                       } from '../../../modules/bwamem2/mem/main'
include { BWA_INDEX                         } from '../../../modules/bwa/index/main.nf'
include { BWA_MEM                           } from '../../../modules/bwa/mem/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../modules/fgbio/callmolecularconsensusreads/main.nf'
include { FGBIO_FASTQTOBAM                  } from '../../../modules/fgbio/fastqtobam/main'
include { FGBIO_GROUPREADSBYUMI             } from '../../../modules/fgbio/groupreadsbyumi/main'
include { SAMBLASTER                        } from '../../../modules/samblaster/main'
include { SAMTOOLS_BAM2FQ                   } from '../../../modules/samtools/bam2fq/main.nf'

workflow CREATE_UMI_CONSENSUS {
    take:
    reads                    // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta                    // channel: [mandatory] /path/to/reference/fasta
    read_structure           // string:  [mandatory] "read_structure"
    groupreadsbyumi_strategy // string:  [mandatory] grouping strategy - default: "Adjacency"
    aligner                  // string:  [mandatory] "bwa-mem" or "bwa-mem2"

    main:
    ch_versions = Channel.empty()

    // Using information in val(read_structure) FASTQ reads are converted into
    // a tagged unmapped BAM file (uBAM)
    FGBIO_FASTQTOBAM ( reads, read_structure )
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.version.first())

    // In order to map uBAM using BWA MEM, we need to convert uBAM to FASTQ
    // but keep the appropriate UMI tags in the FASTQ comment field and produce
    // an interleaved FASTQ file (hence, split = false)
    split = false
    SAMTOOLS_BAM2FQ ( FGBIO_FASTQTOBAM.out.umibam, split )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    // The user can choose here to use either bwa-mem (default) or bwa-mem2
    aligned_bam = Channel.empty()
    if (aligner == "bwa-mem") {
        // Index reference
        BWA_INDEX ( fasta )
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

        // Appropriately tagged, interleaved FASTQ reads are mapped to the reference
        BWA_MEM ( SAMTOOLS_BAM2FQ.out.reads, BWA_INDEX.out.index, false )
        aligned_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    } else {
        // Index reference
        BWAMEM2_INDEX ( fasta )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

        // Appropriately tagged, interleaved FASTQ reads are mapped to the reference
        BWAMEM2_MEM ( SAMTOOLS_BAM2FQ.out.reads, BWAMEM2_INDEX.out.index, false )
        aligned_bam = BWAMEM2_MEM.out.bam
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    }

    // Samblaster is used in order to tag mates information in the BAM file
    // this is used in order to group reads by UMI
    SAMBLASTER ( aligned_bam )
    ch_versions = ch_versions.mix(SAMBLASTER.out.versions.first())

    // Appropriately tagged reads are now grouped by UMI information
    FGBIO_GROUPREADSBYUMI ( SAMBLASTER.out.bam, groupreadsbyumi_strategy )
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())

    // Using the above created groups, a consensus across reads in the same group
    // can be called which will emit a consensus BAM file
    FGBIO_CALLMOLECULARCONSENSUSREADS ( FGBIO_GROUPREADSBYUMI.out.bam )
    ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions.first())

    emit:
    bam_umi       = FGBIO_FASTQTOBAM.out.umibam               // channel: [ val(meta), [ bam ] ]
    bam_umi_group = FGBIO_GROUPREADSBYUMI.out.bam             // channel: [ val(meta), [ bam ] ]
    bam_consensus = FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam // channel: [ val(meta), [ bam ] ]
    versions      = ch_versions                               // channel: [ versions.yml ]
}


include { TABIX_BGZIP as UNZIP_VCFS                     } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_MERGED                   } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'
include { JASMINESV                                     } from '../../../modules/nf-core/jasminesv/main'


workflow VCF_MERGE_JASMINE {
    take:
        ch_vcfs                 // channel: [mandatory] [ meta, vcf ] => The gzipped called VCFs
        ch_bams                 // channel: [optional]  [ meta, bam ] => The BAM files that were used to create the VCFs
        ch_sample_dists         // channel: [optional]  [ meta, txt ] => txt files that contain the distance thresholds for each sample
        ch_fasta                // channel: [optional]  [ fasta ] => The fasta reference file
        ch_fasta_fai            // channel: [optional]  [ fasta_fai ] => The index of the fasta reference file
        ch_chrom_norm           // channel: [optional]  [ txt ] => A txt file containing the chromosomes and their aliases used for normalization
        callers_amount          // value:   [mandatory] callers_amount => The amount of different callers used for the VCFs that need to be merged
        common_meta             // value:   [mandatory] common_meta => A string value stating which meta field contains a common identifier for each sample

    main:

    ch_versions = Channel.empty()

    UNZIP_VCFS(
        ch_vcfs
    )

    ch_versions = ch_versions.mix(UNZIP_VCFS.out.versions.first())

    ch_jasmine_input = UNZIP_VCFS.out.output
        .map { meta, vcf ->
            meta = meta + [id:meta[common_meta]]
            [ meta, vcf ]
        }
        .groupTuple(size:callers_amount)
        .join(ch_bams ?: Channel.value([[],[]]), remainder:true)
        .join(ch_sample_dists ?: Channel.value([[],[]]), remainder:true).view()
        .map { meta, vcfs, bam=null, sample_dists=null ->
            [ meta, vcfs, bam ?: [], sample_dists ?: [] ]
        }
        .filter { meta, vcfs, bam, sample_dists ->
            vcfs != null // Filter out wrongly formed channels after join with remainder true (this won't remove any wanted channels)
        }

    JASMINESV(
        ch_jasmine_input,
        ch_fasta,
        ch_fasta_fai,
        ch_chrom_norm
    )

    ch_versions = ch_versions.mix(JASMINESV.out.versions.first())

    BGZIP_MERGED(
        JASMINESV.out.vcf
    )

    ch_versions = ch_versions.mix(BGZIP_MERGED.out.versions.first())

    TABIX_TABIX(
        BGZIP_MERGED.out.output
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    emit:
    vcf      = BGZIP_MERGED.out.output          // channel: [ val(meta), [ bam ] ]
    tbi      = TABIX_TABIX.out.tbi              // channel: [ val(meta), [ bai ] ]

    versions = ch_versions                      // channel: [ versions.yml ]
}


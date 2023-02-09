include { PREPROCESSGENEANNOTATION  } from '../../../modules/nf-core/aberrantexpression/counting/preprocessgeneannotation'
include { MERGECOUNTS               } from '../../../modules/nf-core/aberrantexpression/counting/mergecounts'
include { FILTERCOUNT               } from '../../../modules/nf-core/aberrantexpression/counting/filtercount'
include { COUNTREADS                } from '../../../modules/nf-core/aberrantexpression/counting/countreads'

include { RUNOUTRIDER               } from '../../../modules/nf-core/aberrantexpression/outrider/runOutrider'
include { RESULTS                   } from '../../../modules/nf-core/aberrantexpression/outrider/outriderresults'

workflow ABERRANT_EXPRESSION {
    take:
        params

    main:
        ch_versions = Channel.empty()


        Channel.from(params.gtf.split(", "))
            .map{iM -> iM.split(":")}
            .set{addProc}
        groups = params.aberrantexpression.groups.split(", ")
        sampleAnn = params.sample_annotation_relative

        // Counting part
        gtfs = addProc.map{iM -> iM[1]}
        anno = addProc.map{iM -> iM[0]}
        PREPROCESSGENEANNOTATION(gtfs)
            .count_reads.map {
                a -> [pathDB: a[0], pathTSV: a[1], count_ranges: a[2]]
            }
            .set {count_reads}

        Channel.fromPath(params.bams).set{ch_bam}
        count_ranges = count_reads.map {it.count_ranges}
        // TODO-csandu: update to reduce
        COUNTREADS(count_ranges, ch_bam)

        // TODO-csandu: Fix dataset
        print(groups)
        MERGECOUNTS(groups, COUNTREADS.out.counts.collect(), count_ranges, sampleAnn)

        txtDB = count_reads.map {it.pathDB}
        FILTERCOUNT(MERGECOUNTS.out.total_counts, txtDB)

        // Outrider part
        RUNOUTRIDER(FILTERCOUNT.out.ods_unfitted,
            params.aberrantexpression.implementation,
            params.aberrantexpression.maxTestedDimensionProportion)

        tsvs = count_reads.map {it.pathTSV}
        RESULTS(params.aberrantexpression.HPO_helper,
            RUNOUTRIDER.out.odsOut,
            tsvs,
            params.aberrantexpression.padjCutoff,
            params.aberrantexpression.zScoreCutoff)
}

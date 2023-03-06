include { PREPROCESSGENEANNOTATION  } from '../../../modules/nf-core/aberrantexpression/counting/preprocessgeneannotation'
include { COUNTREADS                } from '../../../modules/nf-core/aberrantexpression/counting/countreads'
include { MERGECOUNTS               } from '../../../modules/nf-core/aberrantexpression/counting/mergecounts'
include { FILTERCOUNT               } from '../../../modules/nf-core/aberrantexpression/counting/filtercount'

include { RUNOUTRIDER               } from '../../../modules/nf-core/aberrantexpression/outrider/runOutrider'
include { RESULTS                   } from '../../../modules/nf-core/aberrantexpression/outrider/outriderresults'

workflow ABERRANT_EXPRESSION {
    take:
        params

    emit:
        resultsAll

    main:
        ch_versions = Channel.empty()

        PREPROCESSGENEANNOTATION(params.gtfData).result
            .map {it -> ["annotation": it[0], "txtDb": it[1], "geneMap": it[2], "countRanges": it[3]]}
            .set {preprocess_data}

        bams = params.procAnnotation
            .filter (it -> it.RNA_BAM_FILE != "")
        COUNTREADS(preprocess_data, bams.collect()).result
            .map {it -> ["preprocess": it[0], "groups": it[1], "rds": it[2]]}
            .set {countreads_data}

        groups = countreads_data.combine(Channel.from(params.aberrantexpression.groups))
            .filter (it -> it[0].groups.contains(it[1]))
            .map {it -> [it[0].preprocess, it[1], it[0].rds]}
            .groupTuple(by: [0, 1])
            .map {it -> ["preprocess": it[0], "group": it[1], "rds": it[2]]}

        MERGECOUNTS(groups, params.annotation).result
            .map {it -> ["preprocess": it[0], "totalCounts": it[1]]}
            .set {mergecounts_data}

        FILTERCOUNT(mergecounts_data).result
            .map {it -> ["preprocess": it[0], "odsUnfit": it[1]]}
            .set {filter_data}

        // Outrider part
        RUNOUTRIDER(filter_data,
            params.aberrantexpression.implementation,
            params.aberrantexpression.maxTestedDimensionProportion).result
            .map {it -> ["preprocess": it[0], "ods": it[1]]}
            .set {runtoutrider_data}

        RESULTS(
            runtoutrider_data,
            params.hpoFile,
            params.aberrantexpression.padjCutoff,
            params.aberrantexpression.zScoreCutoff)

        resultsAll = RESULTS.out.results_all
}

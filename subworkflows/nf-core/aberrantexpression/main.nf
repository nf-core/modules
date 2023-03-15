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

        /* -------------
         * Counting part
         */
        PREPROCESSGENEANNOTATION(params.gtfData).result
            .map {it -> ["annotation": it[0], "txtDb": it[1], "geneMap": it[2], "countRanges": it[3]]}
            .set {preprocess_data}

        bams = params.procAnnotation
            .filter (it -> it.GENE_COUNTS_FILE == "" && it.SPLICE_COUNTS_DIR == "")
        COUNTREADS(preprocess_data, bams.collect()).result
            .map {it -> ["preprocess": it[0], "groups": it[1], "rds": it[2]]}
            .set {countreads_data}


        /* Concat the exernalData (external_outrider group) to the resulted countreads.
         * The resulted data will contain the RNA_IDs only for external group.
         */
        externalData = params.procAnnotation.filter (it -> it.GENE_COUNTS_FILE != "").combine(preprocess_data)
            .filter {it -> it[0].GENE_ANNOTATION == it[1].annotation}.map {it -> [it[1], it[0].DROP_GROUP, it[0].GENE_COUNTS_FILE, it[0].RNA_ID]}
        groups = countreads_data.combine(Channel.from(params.aberrantexpression.groups))
            .filter (it -> it[0].groups.contains(it[1]))
            .map {it -> [it[0].preprocess, it[1], it[0].rds, "NA"]}
            .concat (externalData)
            .groupTuple(by: [0, 1])
            .map {it -> ["preprocess": it[0], "group": it[1], "rds": it[2], "rnaID": it[3]]}
        MERGECOUNTS(groups, params.annotation).result
            .map {it -> ["preprocess": it[0], "groupName": it[1], "totalCounts": it[2]]}
            .set {mergecounts_data}

        FILTERCOUNT(mergecounts_data, params.aberrantexpression).result
            .map {it -> ["preprocess": it[0], "groupName": it[1], "odsUnfit": it[2]]}
            .set {filter_data}

        /* -------------
         * Outrider part
         */
        RUNOUTRIDER(filter_data,
            params.aberrantexpression.implementation,
            params.aberrantexpression.maxTestedDimensionProportion).result
            .map {it -> ["preprocess": it[0], "groupName": it[1], "ods": it[2]]}
            .set {runoutrider_data}

        RESULTS(runoutrider_data,
            params.hpoFile,
            params.aberrantexpression.padjCutoff,
            params.aberrantexpression.zScoreCutoff).result
            .map {it -> ["preprocess": it[0], "groupName": it[1], "resultsall": it[2], "resultstsv": it[3]]}
            .set {results_data}

        // Prepare the output for moving it to rootDir
        resultsAll = preprocess_data
            .map {it -> [it.annotation, 'general', 'preprocess', [it.txtDb, it.geneMap, it.countRanges]]}
            .concat (countreads_data.map    {it -> [it.preprocess.annotation, 'general', 'countreads', [it.rds]]})
            .concat (mergecounts_data.map   {it -> [it.preprocess.annotation, it.groupName, 'mergecounts', [it.totalCounts]]})
            .concat (filter_data.map        {it -> [it.preprocess.annotation, it.groupName, 'filtercount', [it.odsUnfit]]})
            .concat (runoutrider_data.map   {it -> [it.preprocess.annotation, it.groupName, 'outrider',    [it.ods]]})
            .concat (results_data.map       {it -> [it.preprocess.annotation, it.groupName, 'results',     [it.resultsall, it.resultstsv]]})
}

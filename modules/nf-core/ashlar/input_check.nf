include { SAMPLESHEET_CHECK } from './samplesheet_check'

workflow INPUT_CHECK {

    take:
    samplesheet  // file: ./input_sheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv (header:true, sep:',')
        .map { create_image_channel(it) }
        .set { images }

    emit:
    images                                     // channel: [ val(meta), path(image) ]
    versions = SAMPLESHEET_CHECK.out.versions  // channel: [ versions.yml ]
}

def create_image_channel(LinkedHashMap row) {
    // create meta map
    /*
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    meta.strandedness = row.strandedness

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
    */
    return row
}

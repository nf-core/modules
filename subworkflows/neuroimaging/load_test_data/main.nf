include { TESTDATA_SCILPY      } from '../../../modules/nf-scil/testdata/scilpy/main'

import java.nio.file.Files

workflow LOAD_TEST_DATA {

    take:
    ch_archive
    test_data_prefix

    main:

    ch_versions = Channel.empty()
    test_data_path = Files.createTempDirectory("$test_data_prefix")

    TESTDATA_SCILPY( ch_archive, test_data_path )
    ch_versions = ch_versions.mix(TESTDATA_SCILPY.out.versions.first())

    emit:
    test_data_directory = TESTDATA_SCILPY.out.test_data_directory  // channel: [ test_data_directory ]
    versions            = ch_versions                              // channel: [ versions.yml ]
}

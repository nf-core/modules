#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../../../../modules/nf-core/custom/sratoolsncbisettings/main.nf'

workflow test_sratoolsncbisettings_with_good_existing {

    file(params.settings_path).mkdirs()
    def settings = file(params.test_data['generic']['config']['ncbi_user_settings'], checkIfExists: true)
    settings.copyTo(params.settings_file)

    CUSTOM_SRATOOLSNCBISETTINGS()
}

workflow test_sratoolsncbisettings_with_bad_existing {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = '''
    ## auto-generated configuration file - DO NOT EDIT ##

    config/default = "false"
    /repository/remote/main/CGI/resolver-cgi = "https://trace.ncbi.nlm.nih.gov/Traces/names/names.fcgi"
    /repository/remote/protected/CGI/resolver-cgi = "https://trace.ncbi.nlm.nih.gov/Traces/names/names.fcgi"
    /repository/user/ad/public/apps/file/volumes/flatAd = "."
    /repository/user/ad/public/apps/refseq/volumes/refseqAd = "."
    /repository/user/ad/public/apps/sra/volumes/sraAd = "."
    /repository/user/ad/public/apps/sraPileup/volumes/ad = "."
    /repository/user/ad/public/apps/sraRealign/volumes/ad = "."
    /repository/user/ad/public/apps/wgs/volumes/wgsAd = "."
    /repository/user/ad/public/root = "."
    /repository/user/default-path = "/root/ncbi"
    '''.stripIndent()

    CUSTOM_SRATOOLSNCBISETTINGS()
}

workflow test_sratoolsncbisettings_with_nonexisting {
    def settings = file(params.settings_file)
    settings.delete()

    CUSTOM_SRATOOLSNCBISETTINGS()
}

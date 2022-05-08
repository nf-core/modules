#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRATOOLS_NCBISETTINGS } from '../../../../modules/sratools/ncbisettings/main.nf'

workflow test_sratools_ncbisettings_with_good_existing {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = "/LIBS/GUID = \"5b0d4b7d-88c7-4802-98fd-e3afd06feb32\"\n/libs/cloud/report_instance_identity = \"true\"\n"

    SRATOOLS_NCBISETTINGS()
}

workflow test_sratools_ncbisettings_with_bad_existing {

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

    SRATOOLS_NCBISETTINGS()
}

workflow test_sratools_ncbisettings_with_nonexisting {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.delete()

    SRATOOLS_NCBISETTINGS()
}

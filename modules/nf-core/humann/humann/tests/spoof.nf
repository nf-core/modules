process SPOOF_METAPHLAN_VERSION {
    input:
    tuple val(meta), path(profile)

    output:
    tuple val(meta), path("*_mock_profile.txt"), emit: profile

    script:
    """
    cat ${profile} | sed "s|mpa_vJan21_TOY_CHOCOPhlAnSGB_202103|vOct22_CHOCOPhlAnSGB_202403|g" > test_mock_profile.txt
    """
}

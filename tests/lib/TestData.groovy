class TestData {
    // Same default as params.modules_testdata_base_path in tests/config/nf-test.config
    static String modulesBase() {
        System.getenv('NXF_MODULES_TESTDATA_BASE') ?: 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
    }
}

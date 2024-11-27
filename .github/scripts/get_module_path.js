function get_module_names(filter_modules_files) {
    return [
        ...new Set(
            filter_modules_files.map((path) =>
                path
                    .replace("tests/", "")
                    .replace("modules/nf-core/", "")
                    .split("/")
                    .slice(0, 2)
                    .filter((x) => x !== "main.nf" && x !== "tests" && x !== "meta.yml" && x !== "environment.yml")
                    .join("/")
            )
        ),
    ];
}

test_case_1 = [
    "modules/nf-core/umicollapse/tests/main.nf.test",
    "modules/nf-core/umicollapse/tests/main.nf.test.snap",
    "modules/nf-core/umicollapse/tests/nextflow.config",
    "modules/nf-core/umicollapse/tests/nextflow_PE.config",
    "modules/nf-core/umicollapse/tests/nextflow_SE.config",
    "modules/nf-core/umitools/dedup/tests/main.nf.test",
    "modules/nf-core/umitools/dedup/tests/main.nf.test.snap",
];
result_1 = ["umicollapse", "umitools/dedup"];

console.assert(get_module_names(test_case_1) === result_1, "%o", { "Test Case 1": get_module_names(test_case_1) });

test_case_2 = ["modules/nf-core/umicollapse/tests/main.nf.test"];

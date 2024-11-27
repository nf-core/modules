// Test this locally with:
// node .github/scripts/get_module_path.js

// Core functionality that processes the file paths
function get_module_names(filter_modules_files) {
    return [
        ...new Set(
            filter_modules_files.map((path) =>
                path
                    .replace("tests/", "")
                    .replace("modules/nf-core/", "")
                    .split("/")
                    .slice(0, 2)
                    .filter(
                        (x) =>
                            !x.startsWith("main.nf") &&
                            x !== "tests" &&
                            x !== "meta.yml" &&
                            x !== "environment.yml" &&
                            !x.endsWith(".snap") &&
                            !x.endsWith(".config")
                    )
                    .join("/")
            )
        ),
    ];
}

// GitHub Actions entry point
function run({ github, context }) {
    try {
        // Get the files from the context
        const files = JSON.parse(context.payload.inputs?.files || "[]");
        const result = get_module_names(files);

        // Return json result for GitHub Actions output
        return result;
    } catch (error) {
        console.error("Error processing module paths:", error);
        throw error;
    }
}

// Test cases
function runTests() {
    const test_case_1 = [
        "modules/nf-core/umicollapse/tests/main.nf.test",
        "modules/nf-core/umicollapse/tests/main.nf.test.snap",
        "modules/nf-core/umicollapse/tests/nextflow.config",
        "modules/nf-core/umicollapse/tests/nextflow_PE.config",
        "modules/nf-core/umicollapse/tests/nextflow_SE.config",
        "modules/nf-core/umitools/dedup/tests/main.nf.test",
        "modules/nf-core/umitools/dedup/tests/main.nf.test.snap",
        "modules/nf-core/umitools/dedup/main.nf",
    ];
    const result_1 = ["umicollapse", "umitools/dedup"];

    const test_case_2 = [
        "modules/nf-core/mafft/align/environment.yml",
        "modules/nf-core/mafft/align/main.nf",
        "modules/nf-core/mafft/align/meta.yml",
        "modules/nf-core/mafft/align/tests/main.nf.test",
        "modules/nf-core/mafft/align/tests/main.nf.test.snap",
        "modules/nf-core/mafft/guidetree/environment.yml",
        "modules/nf-core/mafft/guidetree/main.nf",
        "modules/nf-core/mafft/guidetree/meta.yml",
        "modules/nf-core/mafft/guidetree/tests/main.nf.test",
        "modules/nf-core/mafft/guidetree/tests/main.nf.test.snap",
        "tests/modules/nf-core/epang/split/main.nf",
    ];
    const result_2 = ["mafft/align", "mafft/guidetree", "epang/split"];

    console.assert(JSON.stringify(get_module_names(test_case_1)) === JSON.stringify(result_1), "Test case 1 failed");
    console.assert(JSON.stringify(get_module_names(test_case_2)) === JSON.stringify(result_2), "Test case 2 failed");

    console.log("All tests passed!");
}

// Handle different execution contexts
if (require.main === module) {
    // Script was run directly (node get_module_path.js)
    runTests();
} else {
    // Script was imported/required
    module.exports = run;
}

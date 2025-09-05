# ARM Test Skipping

This document describes the simple system for skipping ARM tests in nf-core modules.

## Overview

Some modules may fail to run tests on ARM64 architecture due to various reasons such as:
- Dependencies not available for ARM64
- Architecture-specific compilation issues
- Upstream package limitations

To handle these cases gracefully, we have implemented a simple system that allows temporarily disabling ARM tests for specific modules until the issues are resolved.

## How It Works

### Skip Marker Files

When a module consistently fails ARM tests, you can create a `.skip-arm` file in the module directory. This file:

1. **Disables ARM tests** for that module in CI/CD pipelines
2. **Can contain a link** to the wave build failure (optional)
3. **Can be easily removed** when the issue is resolved

### Workflow Integration

The GitHub Actions workflow automatically:
- Checks for `.skip-arm` files before running tests on ARM runners
- Skips ARM tests for modules with skip markers
- Continues with AMD64 tests as normal
- Logs which tests are being skipped

## Managing ARM Test Skipping

### Simple File Management

**To skip ARM tests for a module:**
```bash
touch modules/nf-core/problematic-module/tests/.skip-arm
```

**To link to a wave build failure (optional):**
```bash
echo "https://wave.seqera.io/build/12345" > modules/nf-core/problematic-module/tests/.skip-arm
```

**To re-enable ARM tests:**
```bash
rm modules/nf-core/problematic-module/tests/.skip-arm
```

**To list modules with ARM tests disabled:**
```bash
find modules/nf-core -name ".skip-arm"
```

## File Format

The `.skip-arm` file can be:
- **Empty**: Just presence of the file disables ARM tests
- **Contain a URL**: Link to wave build failure or GitHub issue
- **Contain notes**: Brief description of the issue

Examples:
```bash
# Empty file
touch modules/nf-core/fastqc/tests/.skip-arm

# With wave build failure link
echo "https://wave.seqera.io/build/12345" > modules/nf-core/fastqc/tests/.skip-arm

# With GitHub issue link
echo "https://github.com/nf-core/modules/issues/12345" > modules/nf-core/fastqc/tests/.skip-arm

# With brief note
echo "bioconda package unavailable for ARM64" > modules/nf-core/fastqc/tests/.skip-arm
```

## Workflow Impact

### Test Filtering

- Modules with `.skip-arm` will have ARM tests skipped
- AMD64 tests continue normally on all profiles (conda, docker, singularity)
- ARM builds still happen when `environment.yml` files are changed (builds are separate from tests)

### CI/CD Logs

When ARM tests are skipped, you'll see logs like:

```
⚠️  Skipping ARM tests for modules/nf-core/fastqc
```

## Best Practices

### When to Add Skip Markers

- ARM tests consistently fail for legitimate architecture reasons
- The failure is due to missing ARM64 packages or architecture-specific issues
- The issue has been investigated and documented
- An upstream issue has been reported (when applicable)

### Documentation Requirements

When adding a skip marker, consider including:
- Link to wave build failure or GitHub issue in the file
- Brief note about the reason (if not obvious from linked issue)

### Regular Review

Periodically review ARM test skips:
- Check if upstream issues have been resolved
- Test if packages have become available for ARM64
- Remove markers when issues are fixed

### Re-enabling ARM Tests

To re-enable ARM tests for a module:
1. Remove the `.skip-arm` file: `rm modules/nf-core/module/.skip-arm`
2. Test the ARM build manually if possible
3. Monitor the CI for successful ARM tests

## Examples

### Example 1: Package Not Available

```bash
# Skip ARM tests due to missing bioconda package
echo "https://github.com/bioconda/bioconda-recipes/issues/12345" > modules/nf-core/tool/tests/.skip-arm
```

### Example 2: Wave Build Failure

```bash
# Skip ARM tests due to wave container build failure
echo "https://wave.seqera.io/build/abc123def456" > modules/nf-core/tool/tests/.skip-arm
```

### Example 3: Simple Skip

```bash
# Just skip without detailed tracking
touch modules/nf-core/tool/tests/.skip-arm
```

## Troubleshooting

### Workflow Still Running ARM Tests

If ARM tests still run despite having a skip marker:
1. Check the file name is exactly `.skip-arm`
2. Verify the file is in the correct module tests directory: `modules/nf-core/modulename/tests/.skip-arm`
3. Check the workflow logs for any errors in the filtering step

### Re-enabling Tests

To test if ARM tests can be re-enabled:
1. Remove the skip marker file: `rm modules/nf-core/module/tests/.skip-arm`
2. Create a small PR that modifies the module
3. Monitor the CI workflow for successful ARM tests
4. If still failing, re-add the marker with updated information

### Finding All Skipped Modules

```bash
# List all modules with ARM tests disabled
find modules/nf-core -name ".skip-arm" -exec dirname {} \;

# See what's in the skip files
find modules/nf-core -name ".skip-arm" -exec echo "=== {} ===" \; -exec cat {} \;
```

# ARM Build Failure Management

This document describes the system for managing ARM build failures in nf-core modules.

## Overview

Some modules may fail to build on ARM64 architecture due to various reasons such as:
- Dependencies not available for ARM64
- Architecture-specific compilation issues
- Upstream package limitations

To handle these cases gracefully, we have implemented a system that allows temporarily disabling ARM builds for specific modules until the issues are resolved.

## How It Works

### Failure Marker Files

When a module consistently fails ARM builds, you can create a `.arm-build-failure.yml` file in the module directory. This file:

1. **Disables ARM builds** for that module in CI/CD pipelines
2. **Documents the failure reason** and relevant information
3. **Can be easily removed** when the issue is resolved

### Workflow Integration

The GitHub Actions workflow automatically:
- Checks for `.arm-build-failure.yml` files before building containers
- Skips ARM builds (`linux/arm64`) for modules with failure markers
- Continues with AMD64 builds (`linux/amd64`) as normal
- Logs the reason for skipping ARM builds

## Managing ARM Build Failures

### Using the Helper Script

We provide a helper script to manage ARM build failure markers:

```bash
# Add a failure marker
./scripts/manage-arm-failures.py add fastqc \
  --reason "Package xyz not available for ARM64" \
  --reported-by "your-github-username" \
  --issue-url "https://github.com/nf-core/modules/issues/12345"

# Check if a module has ARM builds disabled
./scripts/manage-arm-failures.py check fastqc

# List all modules with ARM build failures
./scripts/manage-arm-failures.py list

# Remove a failure marker (re-enable ARM builds)
./scripts/manage-arm-failures.py remove fastqc
```

### Manual Management

You can also manually create/edit the `.arm-build-failure.yml` files:

```yaml
# ARM Build Failure Marker
# This file indicates that ARM builds are disabled for this module
# Remove this file to re-enable ARM builds

reason: "Package XYZ not available for ARM64 architecture"
date: "2024-01-15"
reported_by: "github-username"
issue_url: "https://github.com/nf-core/modules/issues/12345"
error_details: |
  The conda package 'some-package' does not have ARM64 builds available.
  This causes the wave container build to fail consistently.
notes: |
  Alternative solutions investigated:
  - Tried using different base image: failed
  - Contacted upstream maintainer: waiting for response
```

## File Format

The `.arm-build-failure.yml` file supports the following fields:

- **`reason`** (required): Brief description of why ARM builds fail
- **`date`** (required): Date when the failure was first reported
- **`reported_by`** (optional): GitHub username of the person reporting the issue
- **`issue_url`** (optional): Link to the related GitHub issue
- **`error_details`** (optional): Detailed error information or logs
- **`notes`** (optional): Additional notes, investigation details, or workarounds attempted

## Workflow Impact

### Conda-based Builds (`conda-wave` job)

- Modules with `.arm-build-failure.yml` will skip ARM builds
- Docker and Singularity profiles are both affected
- AMD64 builds continue normally

### Dockerfile-based Builds (`dockerfile-wave` job)

- Modules with `.arm-build-failure.yml` will skip ARM builds
- AMD64 builds continue normally

### CI/CD Logs

When ARM builds are skipped, you'll see logs like:

```
⚠️  ARM builds disabled for modules/nf-core/fastqc (failure marker found)
   Reason: Package XYZ not available for ARM64 architecture
   Date: 2024-01-15
   Skipping linux/arm64 build for docker
   Skipping linux/arm64 build for singularity
```

## Best Practices

### When to Add Failure Markers

- ARM builds consistently fail for legitimate architecture reasons
- The failure is not due to temporary issues (network, quota, etc.)
- The issue has been investigated and documented
- An upstream issue has been reported (when applicable)

### Documentation Requirements

When adding a failure marker, please include:
- A clear reason for the failure
- Link to any related GitHub issues
- Details of investigation or workarounds attempted
- Date when the failure was first observed

### Regular Review

Periodically review ARM build failures:
- Check if upstream issues have been resolved
- Test if packages have become available for ARM64
- Remove markers when issues are fixed

### Re-enabling ARM Builds

To re-enable ARM builds for a module:
1. Remove the `.arm-build-failure.yml` file
2. Test the ARM build manually (if possible)
3. Monitor the CI for successful ARM builds

## Examples

### Example 1: Package Not Available

```yaml
reason: "bioconda package 'tool-xyz' has no ARM64 builds"
date: "2024-01-15"
reported_by: "maintainer-username"
issue_url: "https://github.com/bioconda/bioconda-recipes/issues/12345"
error_details: |
  Wave build fails with:
  PackagesNotFoundError: The following packages are not available from current channels:
    - tool-xyz[version='>=1.0.0',build=*linux-aarch64]
notes: |
  Contacted bioconda maintainers about ARM64 support.
  Alternative tools investigated but none suitable.
```

### Example 2: Architecture-Specific Compilation Issue

```yaml
reason: "Compilation fails on ARM64 due to assembly code"
date: "2024-02-01"
reported_by: "developer-username"
issue_url: "https://github.com/upstream/tool/issues/456"
error_details: |
  Build fails during compilation with:
  Error: unsupported instruction on ARM64 architecture
notes: |
  Upstream aware of issue, fix planned for next major release.
  No workaround available currently.
```

## Troubleshooting

### Script Not Working

If the management script doesn't work:
1. Ensure you're in the repository root
2. Check that `uv` is installed
3. Verify the modules directory structure

### Workflow Still Running ARM Builds

If ARM builds still run despite having a failure marker:
1. Check the file name is exactly `.arm-build-failure.yml`
2. Verify the file is in the correct module directory
3. Ensure the YAML syntax is valid
4. Check the workflow logs for any errors in the filtering step

### Re-enabling Builds

To test if ARM builds can be re-enabled:
1. Remove the failure marker file
2. Create a small PR that touches the module's `environment.yml`
3. Monitor the CI workflow for successful ARM builds
4. If still failing, re-add the marker with updated information

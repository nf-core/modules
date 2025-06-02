#!/usr/bin/env uv run
# /// script
# dependencies = [
#   "pyyaml",
# ]
# ///

"""
Script to manage ARM build failure markers.
This script can add, remove, and list ARM failure markers for modules.
"""

import yaml
import os
import sys
import glob
from datetime import datetime
from pathlib import Path


def find_modules():
    """Find all modules and subworkflows."""
    modules = glob.glob("modules/nf-core/*")
    subworkflows = glob.glob("subworkflows/nf-core/*")
    return sorted([m for m in modules + subworkflows if os.path.isdir(m)])


def get_failure_marker_path(module_dir):
    """Get the path to the ARM failure marker for a module."""
    return os.path.join(module_dir, ".arm-build-failure.yml")


def has_failure_marker(module_dir):
    """Check if a module has an ARM failure marker."""
    return os.path.exists(get_failure_marker_path(module_dir))


def create_failure_marker(module_dir, reason="Manual ARM build disable", reported_by="manual"):
    """Create an ARM failure marker for a module."""
    failure_file = get_failure_marker_path(module_dir)

    failure_data = {
        "reason": reason,
        "date": datetime.now().strftime("%Y-%m-%d"),
        "reported_by": reported_by
    }

    with open(failure_file, 'w') as f:
        f.write("# ARM Build Failure Marker\n")
        f.write("# This file indicates that ARM builds are disabled for this module\n")
        f.write("# Remove this file to re-enable ARM builds\n\n")
        yaml.dump(failure_data, f, default_flow_style=False)

    print(f"✅ Created ARM failure marker for {module_dir}")
    return failure_file


def remove_failure_marker(module_dir):
    """Remove an ARM failure marker for a module."""
    failure_file = get_failure_marker_path(module_dir)

    if os.path.exists(failure_file):
        os.remove(failure_file)
        print(f"✅ Removed ARM failure marker for {module_dir}")
        return True
    else:
        print(f"❌ No ARM failure marker found for {module_dir}")
        return False


def list_failure_markers():
    """List all modules with ARM failure markers."""
    modules = find_modules()
    failed_modules = []

    for module_dir in modules:
        if has_failure_marker(module_dir):
            failure_file = get_failure_marker_path(module_dir)
            try:
                with open(failure_file, 'r') as f:
                    failure_info = yaml.safe_load(f) or {}
                failed_modules.append((module_dir, failure_info))
            except Exception as e:
                failed_modules.append((module_dir, {"error": str(e)}))

    if failed_modules:
        print(f"📋 Found {len(failed_modules)} module(s) with ARM failure markers:")
        print()
        for module_dir, failure_info in failed_modules:
            print(f"🔴 {module_dir}")
            if "error" in failure_info:
                print(f"   Error reading marker: {failure_info['error']}")
            else:
                print(f"   Reason: {failure_info.get('reason', 'No reason specified')}")
                print(f"   Date: {failure_info.get('date', 'No date specified')}")
                print(f"   Reported by: {failure_info.get('reported_by', 'Unknown')}")
            print()
    else:
        print("✅ No modules with ARM failure markers found")

    return failed_modules


def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) < 2:
        print("Usage: manage_arm_failures.py <command> [args]")
        print("\nCommands:")
        print("  list                     - List all modules with ARM failure markers")
        print("  add <module_dir> [reason] - Add ARM failure marker to module")
        print("  remove <module_dir>      - Remove ARM failure marker from module")
        print("  clean                    - Remove all ARM failure markers")
        print("\nExamples:")
        print("  manage_arm_failures.py list")
        print("  manage_arm_failures.py add modules/nf-core/fastqc 'Known ARM issue'")
        print("  manage_arm_failures.py remove modules/nf-core/fastqc")
        sys.exit(1)

    command = sys.argv[1]

    if command == "list":
        list_failure_markers()

    elif command == "add":
        if len(sys.argv) < 3:
            print("Error: Module directory required")
            sys.exit(1)

        module_dir = sys.argv[2]
        reason = sys.argv[3] if len(sys.argv) > 3 else "Manual ARM build disable"

        if not os.path.isdir(module_dir):
            print(f"Error: Directory {module_dir} does not exist")
            sys.exit(1)

        create_failure_marker(module_dir, reason, "manual")

    elif command == "remove":
        if len(sys.argv) < 3:
            print("Error: Module directory required")
            sys.exit(1)

        module_dir = sys.argv[2]
        remove_failure_marker(module_dir)

    elif command == "clean":
        modules = find_modules()
        removed_count = 0

        for module_dir in modules:
            if remove_failure_marker(module_dir):
                removed_count += 1

        print(f"✅ Removed {removed_count} ARM failure marker(s)")

    else:
        print(f"Error: Unknown command '{command}'")
        print("Use 'list', 'add', 'remove', or 'clean'")
        sys.exit(1)


if __name__ == "__main__":
    main()

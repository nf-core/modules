#!/usr/bin/env python3
"""Add arm_failure tag to nf-test files for modules that fail on ARM"""

import sys
import re
from pathlib import Path

def add_arm_failure_tag(test_file_path):
    """Add arm_failure tag to the test file if it doesn't already exist"""

    if not test_file_path.exists():
        print(f"Test file not found: {test_file_path}")
        return False

    content = test_file_path.read_text()

    # Check if arm_failure tag already exists
    if 'tag "arm_failure"' in content:
        print(f"ARM failure tag already exists in {test_file_path}")
        return False

    # Find the last tag line in the file - looking for pattern: tag "something"
    tag_pattern = r'(\s*tag\s+"[^"]+"\s*\n)'
    tag_matches = list(re.finditer(tag_pattern, content))

    if tag_matches:
        # Get the last tag match
        last_tag = tag_matches[-1]
        # Insert the new tag after the last existing tag with proper formatting
        before = content[:last_tag.end()]
        after = content[last_tag.end():]
        new_content = before + '    tag "arm_failure"\n\n' + after
    else:
        # No existing tags, add after the process line
        process_pattern = r'(process\s+"[^"]+"\s*\n)'
        process_match = re.search(process_pattern, content)
        if process_match:
            before = content[:process_match.end()]
            after = content[process_match.end():]
            new_content = before + '\n    tag "arm_failure"\n' + after
        else:
            print(f"Could not find process line in {test_file_path}")
            return False

    # Write the updated content
    test_file_path.write_text(new_content)
    print(f"Added ARM failure tag to {test_file_path}")
    return True

def main():
    if len(sys.argv) < 2:
        print("Usage: add_arm_failure_tag.py <module_path> [module_path2] ...")
        print("Example: add_arm_failure_tag.py modules/nf-core/fastqc")
        sys.exit(1)

    success_count = 0
    for module_path in sys.argv[1:]:
        module_dir = Path(module_path)
        test_file = module_dir / "tests" / "main.nf.test"

        if add_arm_failure_tag(test_file):
            success_count += 1

    print(f"Successfully added ARM failure tags to {success_count} module(s)")

if __name__ == "__main__":
    main()

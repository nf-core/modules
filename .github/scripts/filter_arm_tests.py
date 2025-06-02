#!/usr/bin/env python
"""Skip ARM tests for modules with .skip-arm marker files"""
import json, sys, os

def main():
    if len(sys.argv) != 2:
        print(json.dumps([]))
        return

    paths = json.loads(sys.argv[1])
    filtered = []

    for path in paths:
        # Extract module dir (modules/nf-core/xxx or subworkflows/nf-core/xxx)
        parts = path.split("/")
        if len(parts) >= 3 and parts[1] == "nf-core":
            module_dir = "/".join(parts[:3])
            if os.path.exists(f"{module_dir}/.skip-arm"):
                print(f"⚠️  Skipping ARM tests for {path}", file=sys.stderr)
                continue
        filtered.append(path)

    print(json.dumps(filtered))

if __name__ == "__main__": main()

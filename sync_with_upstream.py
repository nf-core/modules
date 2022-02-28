#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import tempfile

import tolit.slack

script_description = """
Script to merge the upstream nf-core 'master' branch onto our own nf-core-modules ('main' branch)
"""

GITHUB_SANGER_MODULES = 'git@github.com:sanger-tol/nf-core-modules.git'
GITHUB_NFCORE_MODULES = 'git@github.com:nf-core/modules.git'
SANGER_REMOTE_NAME = 'sanger'
NFCORE_REMOTE_NAME = 'upstream'
SANGER_BRANCH_NAME = 'main'
NFCORE_BRANCH_NAME = 'master'


def run_and_check_command(args):
    print("Running", args)
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode:
        msg_items = [f'The command `{" ".join(args)}` failed with return code {proc.returncode}']
        if proc.stdout:
            msg_items.extend([
                '_stdout_:',
                '```',
                proc.stdout,
                '```',
            ])
        if proc.stderr:
            msg_items.extend([
                '_stderr_:',
                '```',
                proc.stderr,
                '```',
            ])
        tolit.slack.report_to_slack(os.environ["SLACK_WEBHOOK_PERSO"], "\n".join(msg_items))
        sys.exit(1)
    else:
        print("OK")
        print(proc.stdout)

def main():
    parser = argparse.ArgumentParser(description=script_description)
    parser.add_argument('-n', '--dry-run', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmpdirname:
        run_and_check_command(["git", "clone", "-o", SANGER_REMOTE_NAME, "-b", SANGER_BRANCH_NAME, GITHUB_SANGER_MODULES, tmpdirname])
        os.chdir(tmpdirname)
        run_and_check_command(["git", "remote", "add", NFCORE_REMOTE_NAME, GITHUB_NFCORE_MODULES])
        run_and_check_command(["git", "fetch", NFCORE_REMOTE_NAME])
        run_and_check_command(["git", "merge", "--ff", "--no-edit", f"{NFCORE_REMOTE_NAME}/{NFCORE_BRANCH_NAME}"])
        if args.dry_run:
            run_and_check_command(["git", "push", "-n", SANGER_REMOTE_NAME])
        else:
            run_and_check_command(["git", "push", SANGER_REMOTE_NAME])

if __name__ == "__main__":
    main()

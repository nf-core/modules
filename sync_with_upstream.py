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

# Define a cron job such as:
#30 2 * * * /bin/bash -cl 'env PYTHONPATH=$HOME/workspace/tol-it/lib /software/python-3.9.2/bin/python3 $HOME/workspace/tol-it/nextflow/modules/sync_with_upstream.py > workspace/tol-it/nextflow/modules/sync_with_upstream.cron.txt 2>&1'


GITHUB_SANGER_MODULES = 'git@github.com:sanger-tol/nf-core-modules.git'
GITHUB_NFCORE_MODULES = 'git@github.com:nf-core/modules.git'
SANGER_REMOTE_NAME = 'sanger'
NFCORE_REMOTE_NAME = 'upstream'
SANGER_BRANCH_NAME = 'main'
NFCORE_BRANCH_NAME = 'master'


# Run the command, and report failures to Slack
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
    # No parameters are exposed. Only a dry-run option for testing
    parser = argparse.ArgumentParser(description=script_description)
    parser.add_argument('-n', '--dry-run', action='store_true')
    args = parser.parse_args()
    # Work in an isolated directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        # A bunch of commands we run one after the other
        run_and_check_command(["git", "clone", "-o", SANGER_REMOTE_NAME, "-b", SANGER_BRANCH_NAME, GITHUB_SANGER_MODULES, tmpdirname])
        os.chdir(tmpdirname)
        run_and_check_command(["git", "remote", "add", NFCORE_REMOTE_NAME, GITHUB_NFCORE_MODULES])
        run_and_check_command(["git", "fetch", NFCORE_REMOTE_NAME])
        run_and_check_command(["git", "merge", "--no-edit", f"{NFCORE_REMOTE_NAME}/{NFCORE_BRANCH_NAME}"])
        if args.dry_run:
            run_and_check_command(["git", "push", "-n", SANGER_REMOTE_NAME])
            run_and_check_command(["git", "push", "-n", SANGER_REMOTE_NAME, f"{NFCORE_REMOTE_NAME}/{NFCORE_BRANCH_NAME}:nf-core_master"])
        else:
            run_and_check_command(["git", "push", SANGER_REMOTE_NAME])
            run_and_check_command(["git", "push", SANGER_REMOTE_NAME, f"{NFCORE_REMOTE_NAME}/{NFCORE_BRANCH_NAME}:nf-core_master"])

if __name__ == "__main__":
    main()

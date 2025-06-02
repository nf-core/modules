#!/usr/bin/env uv run
# /// script
# dependencies = [
#   "pyyaml",
#   "click",
# ]
# ///

"""
Helper script to manage ARM build failure markers for nf-core modules.

This script allows you to:
- Add ARM build failure markers to modules
- Remove ARM build failure markers from modules
- List all modules with ARM build failures
- Check if a specific module has ARM build failures
"""

import click
import yaml
import os
from pathlib import Path
from datetime import datetime


FAILURE_MARKER_FILENAME = ".arm-build-failure.yml"


def get_modules_dir():
    """Get the modules directory path."""
    script_dir = Path(__file__).parent
    modules_dir = script_dir.parent / "modules" / "nf-core"
    if not modules_dir.exists():
        raise click.ClickException(f"Modules directory not found: {modules_dir}")
    return modules_dir


def get_module_path(module_name):
    """Get the path to a specific module."""
    modules_dir = get_modules_dir()
    module_path = modules_dir / module_name
    if not module_path.exists():
        raise click.ClickException(f"Module not found: {module_name}")
    return module_path


def get_failure_marker_path(module_name):
    """Get the path to the ARM failure marker file for a module."""
    module_path = get_module_path(module_name)
    return module_path / FAILURE_MARKER_FILENAME


@click.group()
def cli():
    """Manage ARM build failure markers for nf-core modules."""
    pass


@cli.command()
@click.argument('module_name')
@click.option('--reason', '-r', required=True, help='Reason for ARM build failure')
@click.option('--reported-by', '-u', help='GitHub username of reporter')
@click.option('--issue-url', '-i', help='URL to related GitHub issue')
@click.option('--error-details', '-e', help='Detailed error information')
@click.option('--notes', '-n', help='Additional notes')
def add(module_name, reason, reported_by, issue_url, error_details, notes):
    """Add an ARM build failure marker to a module."""
    failure_marker_path = get_failure_marker_path(module_name)

    if failure_marker_path.exists():
        if not click.confirm(f"ARM failure marker already exists for {module_name}. Overwrite?"):
            return

    failure_data = {
        'reason': reason,
        'date': datetime.now().strftime('%Y-%m-%d'),
    }

    if reported_by:
        failure_data['reported_by'] = reported_by
    if issue_url:
        failure_data['issue_url'] = issue_url
    if error_details:
        failure_data['error_details'] = error_details
    if notes:
        failure_data['notes'] = notes

    with open(failure_marker_path, 'w') as f:
        f.write("# ARM Build Failure Marker\n")
        f.write("# This file indicates that ARM builds are disabled for this module\n")
        f.write("# Remove this file to re-enable ARM builds\n\n")
        yaml.dump(failure_data, f, default_flow_style=False, sort_keys=False)

    click.echo(f"✅ Added ARM failure marker for {module_name}")


@cli.command()
@click.argument('module_name')
def remove(module_name):
    """Remove an ARM build failure marker from a module."""
    failure_marker_path = get_failure_marker_path(module_name)

    if not failure_marker_path.exists():
        click.echo(f"❌ No ARM failure marker found for {module_name}")
        return

    if click.confirm(f"Remove ARM failure marker for {module_name}?"):
        failure_marker_path.unlink()
        click.echo(f"✅ Removed ARM failure marker for {module_name}")


@cli.command()
@click.argument('module_name')
def check(module_name):
    """Check if a module has an ARM build failure marker."""
    failure_marker_path = get_failure_marker_path(module_name)

    if failure_marker_path.exists():
        click.echo(f"⚠️  {module_name} has ARM builds disabled")

        try:
            with open(failure_marker_path, 'r') as f:
                failure_data = yaml.safe_load(f)
                if failure_data:
                    click.echo(f"   Reason: {failure_data.get('reason', 'Not specified')}")
                    click.echo(f"   Date: {failure_data.get('date', 'Not specified')}")
                    if 'reported_by' in failure_data:
                        click.echo(f"   Reported by: {failure_data['reported_by']}")
                    if 'issue_url' in failure_data:
                        click.echo(f"   Issue: {failure_data['issue_url']}")
        except Exception as e:
            click.echo(f"   Error reading failure data: {e}")
    else:
        click.echo(f"✅ {module_name} has ARM builds enabled")


@cli.command()
def list():
    """List all modules with ARM build failure markers."""
    modules_dir = get_modules_dir()
    failed_modules = []

    for module_dir in modules_dir.iterdir():
        if module_dir.is_dir():
            failure_marker = module_dir / FAILURE_MARKER_FILENAME
            if failure_marker.exists():
                try:
                    with open(failure_marker, 'r') as f:
                        failure_data = yaml.safe_load(f)
                        failed_modules.append({
                            'name': module_dir.name,
                            'reason': failure_data.get('reason', 'Not specified') if failure_data else 'Not specified',
                            'date': failure_data.get('date', 'Not specified') if failure_data else 'Not specified'
                        })
                except Exception as e:
                    failed_modules.append({
                        'name': module_dir.name,
                        'reason': f'Error reading file: {e}',
                        'date': 'Unknown'
                    })

    if not failed_modules:
        click.echo("✅ No modules have ARM build failures")
        return

    click.echo(f"⚠️  Found {len(failed_modules)} modules with ARM build failures:\n")

    for module in sorted(failed_modules, key=lambda x: x['name']):
        click.echo(f"  {module['name']}")
        click.echo(f"    Reason: {module['reason']}")
        click.echo(f"    Date: {module['date']}")
        click.echo()


if __name__ == '__main__':
    cli()

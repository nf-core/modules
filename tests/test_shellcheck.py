import json
import os
import re
import subprocess
import sys
import textwrap

from tqdm.auto import tqdm

def extract_shell_scripts(nf_content):
    # Match script or shell blocks with triple quotes, excluding groovy in between
    # (e.g. // comments or def variable definitions)
    pattern = re.compile(r'\n\s*(script|shell):\s*(?:.*?\n)*?\s*"""(.*?)"""', re.DOTALL)
    matches = pattern.findall(nf_content)

    scripts = []
    positions = []

    for match in matches:
        script = match[1]
        start_pos = nf_content.find(script)
        scripts.append(script)
        positions.append(start_pos)

    return scripts, positions

def save_scripts(scripts, base_filename):
    script_files = []
    for i, script in enumerate(scripts):
        filename = f"{base_filename}_script_{i+1}.sh"
        with open(filename, mode='w', newline="\n") as f:
            # backslashes are escaped, so we unescape them:
            script = script.replace("\\\\", "\\")
            # the script is indented, we unindent it:
            script = textwrap.dedent(script)
            # Now it looks like the bash script, write it:
            f.write(script)
        script_files.append(filename)
    return script_files

def run_shellcheck(script_files, nf_content, positions, nf_file):
    markdown_output = []
    json_output = []
    code_block = False  # Keep track of code blocks in markdown file
    for i, script_file in enumerate(script_files):
        with open(script_file) as f:
            first_line = f.readline().strip()
            if first_line.startswith('#!') and 'bash' not in first_line:
                markdown_output.append(f"Skipping {script_file} due to non-Bash shebang: {first_line}\n")
                json_output.append({"file": script_file, "message": f"Skipping due to non-Bash shebang: {first_line}"})
                print(f"Skipping {script_file} due to non-Bash shebang: {first_line}")
                continue

        # Let shellcheck ignore SC2154 (ignore using undefined variables). There is a lot of
        # false positives due to variables defined by groovy before the script.
        result = subprocess.run(['shellcheck', '--shell=bash', '--exclude=SC2154', script_file], capture_output=True, text=True)
        start_pos = positions[i]

        if result.stdout or result.stderr:
            if code_block:
                markdown_output.append("\n```\n")
                code_block = False
            markdown_output.append(f"\n## ShellCheck results for `{nf_file}`:\n")
            json_output.append({"file": nf_file, "results": [{"line": -1, "message": []}]})
            print(f"\nShellCheck results for {nf_file}:\n")
            for line in result.stdout.splitlines():
                match = re.match('In (.*\\.sh) line (\\d+):', line)
                if match:
                    script_line_num = int(match.group(2))
                    nf_line_num = nf_content.count('\n', 0, start_pos) + script_line_num
                    if code_block:
                        markdown_output.append("\n```\n")
                        code_block = False
                    markdown_output.append(f"### In `{nf_file}` line {nf_line_num}:\n")
                    markdown_output.append("\n```\n")
                    code_block = True
                    json_output[-1]["results"].append({"line": nf_line_num, "message": [line]})
                    print(f"In {nf_file} line {nf_line_num}:")
                else:
                    markdown_output.append(line + "\n")
                    json_output[-1]["results"][-1]["message"].append(line)
                    print(line)

            if result.stderr:
                markdown_output.append(result.stderr + "\n")
                json_output[-1]["results"].append({"error": result.stderr})
                print(result.stderr)
        os.remove(script_file)

    if code_block:
        markdown_output.append("\n```\n")
        code_block = False
    return (markdown_output, json_output)



def process_nf_file(nf_file):
    with open(nf_file) as f:
        nf_content = f.read().replace("\r\n", "\n")

    scripts, positions = extract_shell_scripts(nf_content)
    base_filename = os.path.splitext(nf_file.replace('/', '_'))[0]
    script_files = save_scripts(scripts, base_filename)
    return run_shellcheck(script_files, nf_content, positions, nf_file)

def main():
    if len(sys.argv) == 1:
        nf_files = []
        for root, dirs, files in os.walk('.'):
            for file in files:
                if file.endswith(".nf"):
                    nf_file = os.path.join(root, file)
                    nf_files.append(nf_file)
    else:
        nf_files = [x for x in sys.argv if x.endswith(".nf")]
    markdown_output = []
    json_output = []
    for i, nf_file in enumerate(tqdm(nf_files, desc="Shellchecking modified .nf files")):
        (md_i, js_i) = process_nf_file(nf_file)
        markdown_output.extend(md_i)
        json_output.extend(js_i)

    if markdown_output:
        with open("shellcheck_output.md", "w") as md_file:
            md_file.writelines(markdown_output)

    if json_output:
        with open("shellcheck_output.json", "w") as json_file:
            json.dump(json_output, json_file, indent=4)

if __name__ == "__main__":
    main()

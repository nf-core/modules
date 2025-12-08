# Copyright 2022 CRS4.
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.

from pathlib import PurePosixPath
from urllib.parse import unquote, urlparse


def get_ci_wf_endpoint(repo_url, ci_wf_name):
    repo_path = PurePosixPath(urlparse(unquote(repo_url)).path)
    if len(repo_path.parts) != 3:  # first one is '/'
        raise ValueError("repository url must be like https://github.com/<OWNER>/<REPO>")
    owner, repo_name = repo_path.parts[-2:]
    return f"repos/{owner}/{repo_name}/actions/workflows/{ci_wf_name}"

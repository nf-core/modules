# Copyright 2019-2025 The University of Manchester, UK
# Copyright 2020-2025 Vlaams Instituut voor Biotechnologie (VIB), BE
# Copyright 2020-2025 Barcelona Supercomputing Center (BSC), ES
# Copyright 2020-2025 Center for Advanced Studies, Research and Development in Sardinia (CRS4), IT
# Copyright 2022-2025 École Polytechnique Fédérale de Lausanne, CH
# Copyright 2024-2025 Data Centre, SciLifeLab, SE
# Copyright 2024-2025 National Institute of Informatics (NII), JP
# Copyright 2025 Senckenberg Society for Nature Research (SGN), DE
# Copyright 2025 European Molecular Biology Laboratory (EMBL), Heidelberg, DE
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from .contextentity import ContextEntity


class TestService(ContextEntity):

    def _empty(self):
        return {
            "@id": self.id,
            "@type": 'TestService'
        }

    @property
    def _default_type(self):
        return "TestService"

    @property
    def name(self):
        return self.get("name")

    @name.setter
    def name(self, name):
        self["name"] = name

    @property
    def url(self):
        return self.get("url")

    @url.setter
    def url(self, url):
        self["url"] = url


JENKINS_ID = "https://w3id.org/ro/terms/test#JenkinsService"
TRAVIS_ID = "https://w3id.org/ro/terms/test#TravisService"
GITHUB_ID = "https://w3id.org/ro/terms/test#GithubService"


def jenkins(crate):
    return TestService(crate, identifier=JENKINS_ID, properties={
        "name": "Jenkins",
        "url": {
            "@id": "https://www.jenkins.io"
        },
    })


def travis(crate):
    return TestService(crate, identifier=TRAVIS_ID, properties={
        "name": "Travis CI",
        "url": {
            "@id": "https://www.travis-ci.com"
        },
    })


def github(crate):
    return TestService(crate, identifier=GITHUB_ID, properties={
        "name": "Github Actions",
        "url": {
            "@id": "https://github.com"
        },
    })


SERVICE_MAP = {
    "jenkins": jenkins,
    "travis": travis,
    "github": github,
}


def get_service(crate, name):
    try:
        func = SERVICE_MAP[name.lower()]
    except KeyError:
        raise ValueError(f"Unknown service: {name}")
    return func(crate)

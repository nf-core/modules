#!/usr/bin/env python

## Copyright 2018-2020 Stian Soiland-Reyes, The University of Manchester, UK
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

"""
Generate arcp URIs with various prefixes.

As detailed in draft-soilandreyes-arcp_, the choice of 
arcp _prefix_ depends on the uniqueness constraints required
to identify the archive.

:func:`arcp_random()` can be used for a fresh arcp URI 
based on a pseudo-random generator. Use 
:func:`urllib.parse.urljoin()` to 
resolve paths within the same archive.

:func:`arcp_uuid()` can be used with a pre-made UUID instance,
for instance loaded from an archive's manifest
or generated with :func:`uuid.uuid4()`

:func:`arcp_location()` can be used to identify an archive based on
its location URL, facilitating a UUID v5 authority.

:func:`arcp_name()` can be used to identify an archive based on its
absolute DNS name or package name within an installation.

.. _draft-soilandreyes-arcp: https://tools.ietf.org/id/draft-soilandreyes-arcp-03.html
"""

__author__      = "Stian Soiland-Reyes <http://orcid.org/0000-0001-9842-9718>"
__copyright__   = "Copyright 2018 The University of Manchester"
__license__     = "Apache License, version 2.0 (https://www.apache.org/licenses/LICENSE-2.0)"

from uuid import uuid4, uuid5, UUID, NAMESPACE_URL

try:
    from urllib.parse import urlunsplit
except:
    from urlparse import urlunsplit

import re
from hashlib import sha256
from base64 import urlsafe_b64encode, urlsafe_b64decode

SCHEME="arcp"

def _reg_name_regex():
    """Compile regular expression for RFC3986_ reg-name production

    _RFC3986: https://www.ietf.org/rfc/rfc3986
    """
    # unreserved    = ALPHA / DIGIT / "-" / "." / "_" / "~"
    unreserved = r"[A-Za-z0-9-._~]"

    # pct-encoded = "%" HEXDIG HEXDIG
    pct_encoded = r"%[0-9A-Fa-f][0-9A-Fa-f]"
    
    # "!" / "$" / "&" / "'" / "(" / ")"
    #  / "*" / "+" / "," / ";" / "="
    sub_delims = r"[!$&'()*+,;=]"

    # reg-name    = *( unreserved / pct-encoded / sub-delims )    
    reg_name = r"^(" + unreserved + r"|" + pct_encoded + sub_delims + r")*$"
    return re.compile(reg_name)
_REG_NAME = _reg_name_regex()

def arcp_uuid(uuid, path="/", query=None, fragment=None):
    """Generate an arcp URI for the given uuid.

    Parameters:
      - uuid -- a uuid string or UUID instance identifying the archive, e.g. ``58ca7fa6-be2f-48e4-8b69-e63fb0d929fe``
      - path -- Optional path within archive.
      - query -- Optional query component.
      - fragment -- Optional fragment component.
    """
    if not isinstance(uuid, UUID):
        # ensure valid UUID
        uuid = UUID(uuid)

    # TODO: Ensure valid path?    
    path = path or ""
    authority = "uuid,%s" % uuid
    s = (SCHEME, authority, path, query, fragment)
    return urlunsplit(s)

def arcp_random(path="/", query=None, fragment=None, uuid=None):
    """Generate an arcp URI using a random uuid.

    Parameters:
      - path -- Optional path within archive.
      - query -- Optional query component.
      - fragment -- Optional fragment component.
      - uuid -- optional UUID v4 string or UUID instance
    """
    if uuid is None:
        uuid = uuid4()
    elif not isinstance(uuid, UUID):
        # ensure valid UUID
        uuid = UUID(uuid)
    if not uuid.version == 4:
        raise Exception("UUID is not v4" % uuid)
    return arcp_uuid(uuid, path=path, query=query, fragment=fragment)

def arcp_location(location, path="/", query=None, fragment=None, namespace=NAMESPACE_URL):
    """Generate an arcp URI for a given archive location.

    Parameters:
      - location: URL or location of archive, e.g. ``http://example.com/data.zip``
      - path -- Optional path within archive.
      - query -- Optional query component.
      - fragment -- Optional fragment component.
      - namespace -- optional namespace UUID for non-URL location.
    """
    # TODO: Ensure location is valid url if NAMESPACE_URL?
    uuid = uuid5(namespace, location)
    return arcp_uuid(uuid, path=path, query=query, fragment=fragment)
    
def arcp_name(name, path="/", query=None, fragment=None):
    """Generate an arcp URI for a given archive name.

    Parameters:
      - name -- Absolute DNS or package name, e.g. ``app.example.com``
      - path -- Optional path within archive.
      - query -- Optional query component.
      - fragment -- Optional fragment component.
      - namespace -- optional namespace UUID for non-URL location.
    """
    if not _REG_NAME.match(name):
        raise Exception("Invalid name: %s" % name)
    authority = "name," + name
    s = (SCHEME, authority, path, query, fragment)
    return urlunsplit(s)

def arcp_hash(bytes=b"", path="/", query=None, fragment=None, hash=None):
    """Generate an arcp URI for a given archive hash checksum.

    Parameters:
      - bytes -- Optional bytes of archive to checksum
      - path -- Optional path within archive.
      - query -- Optional query component.
      - fragment -- Optional fragment component.
      - hash -- Optional hash instance from :func:`hashlib.sha256()`
    
    Either ``bytes`` or ``hash`` must be provided. 
    The ``hash`` parameter can be provided to avoid representing 
    the whole archive bytes in memory.
    """
    if hash is None:
        hash = sha256()
    elif hash.name != "sha256":
        # TODO: Map Python's hash-names to RFC6920
        raise Exception("hash method %s unsupported, try sha256" % hash.name)
    hashmethod = "sha-256"

    # Tip: if bytes == b"" then provided hash param is unchanged
    hash.update(bytes)

    # RFC6920-style hash encoding
    digestB64 = urlsafe_b64encode(hash.digest())
    digestB64 = digestB64.decode("ascii").strip("=")
    authority = "ni,%s;%s" % (hashmethod, digestB64)
    s = (SCHEME, authority, path, query, fragment)
    return urlunsplit(s)


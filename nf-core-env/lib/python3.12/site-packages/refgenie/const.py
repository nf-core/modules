"""
Constant variables for refgenie package. Ones that are integral to refgenconf
and/or refgenieserver should be defined in refgenconf.const
"""
from refgenconf.const import *

PKG_NAME = "refgenie"

BUILD_CMD = "build"
INIT_CMD = "init"
PULL_CMD = "pull"
POPULATE_CMD = "populate"
POPULATE_REMOTE_CMD = "populater"
LIST_LOCAL_CMD = "list"
LIST_REMOTE_CMD = "listr"
GET_ASSET_CMD = "seek"
GET_REMOTE_ASSET_CMD = "seekr"
INSERT_CMD = "add"
REMOVE_CMD = "remove"
GETSEQ_CMD = "getseq"
TAG_CMD = "tag"
ID_CMD = "id"
SUBSCRIBE_CMD = "subscribe"
UNSUBSCRIBE_CMD = "unsubscribe"
ALIAS_CMD = "alias"
COMPARE_CMD = "compare"
UPGRADE_CMD = "upgrade"

GENOME_ONLY_REQUIRED = [REMOVE_CMD, GETSEQ_CMD]

# For each asset we assume a genome is also required
ASSET_REQUIRED = [
    PULL_CMD,
    GET_ASSET_CMD,
    GET_REMOTE_ASSET_CMD,
    BUILD_CMD,
    INSERT_CMD,
    TAG_CMD,
    ID_CMD,
]

SUBPARSER_MESSAGES = {
    INIT_CMD: "Initialize a genome configuration.",
    LIST_LOCAL_CMD: "List available local assets.",
    LIST_REMOTE_CMD: "List available remote assets.",
    PULL_CMD: "Download assets.",
    BUILD_CMD: "Build genome assets.",
    GET_ASSET_CMD: "Get the path to a local asset.",
    GET_REMOTE_ASSET_CMD: "Get the path to a remote asset.",
    INSERT_CMD: "Add local asset to the config file.",
    REMOVE_CMD: "Remove a local asset.",
    GETSEQ_CMD: "Get sequences from a genome.",
    TAG_CMD: "Tag an asset.",
    ID_CMD: "Return the asset digest.",
    SUBSCRIBE_CMD: "Add a refgenieserver URL to the config.",
    UNSUBSCRIBE_CMD: "Remove a refgenieserver URL from the config.",
    ALIAS_CMD: "Interact with aliases.",
    COMPARE_CMD: "Compare two genomes.",
    UPGRADE_CMD: "Upgrade config. This will alter the files on disk.",
    POPULATE_CMD: "Populate registry paths with local paths.",
    POPULATE_REMOTE_CMD: "Populate registry paths with remote paths.",
}


ALIAS_GET_CMD = "get"
ALIAS_SET_CMD = "set"
ALIAS_REMOVE_CMD = "remove"

ALIAS_SUBPARSER_MESSAGES = {
    ALIAS_REMOVE_CMD: "Remove aliases.",
    ALIAS_SET_CMD: "Set aliases.",
    ALIAS_GET_CMD: "Get aliases.",
}

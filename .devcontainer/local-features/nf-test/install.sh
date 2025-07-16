#!/usr/bin/env bash

# download and install
wget -qO- https://get.nf-test.com | bash

# move to location on the path
mv nf-test "$HOME/.local/bin

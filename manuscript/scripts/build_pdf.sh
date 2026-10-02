#!/usr/bin/env bash
# Current Article + SI; historical submission packages are intentionally untouched.
set -euo pipefail
exec bash "$(dirname "$0")/rebuild_current.sh"

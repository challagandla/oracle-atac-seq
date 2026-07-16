#!/usr/bin/env bash
# Backward-compatible entry point. setup.sh is the authoritative installer.
set -Eeuo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
printf '[INFO] scripts/setup_envs.sh now delegates to setup.sh\n'
exec bash "$ROOT_DIR/setup.sh" "$@"

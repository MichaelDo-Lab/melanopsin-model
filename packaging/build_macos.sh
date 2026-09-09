#!/usr/bin/env bash
# Build the Melanopsin Model one-click macOS .app for this machine's architecture.
#
# Creates (or reuses) an isolated virtual environment, installs the pinned
# build dependencies, and runs PyInstaller against packaging/melanopsin_gui.spec.
# The resulting bundle is written to
# dist/MelanopsinModel-v<version>-macos-<arch>.app and zipped beside it, where
# <version> comes from myutils/_version.py and <arch> is arm64 or x86_64.
#
# Usage (from the repository root):
#   bash packaging/build_macos.sh --clean
#   bash packaging/build_macos.sh --skip-venv
set -euo pipefail

CLEAN=0
SKIP_VENV=0
for arg in "$@"; do
    case "$arg" in
        --clean)
            CLEAN=1
            ;;
        --skip-venv)
            SKIP_VENV=1
            ;;
        -h|--help)
            echo "Usage: bash packaging/build_macos.sh [--clean] [--skip-venv]"
            echo "Builds dist/MelanopsinModel-v<version>-macos-<arch>.app for this Mac."
            exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Usage: bash packaging/build_macos.sh [--clean] [--skip-venv]" >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"
echo "Repository root: $REPO_ROOT"

if [[ "$CLEAN" -eq 1 ]]; then
    echo "Cleaning build/ and dist/ ..."
    rm -rf "$REPO_ROOT/build" "$REPO_ROOT/dist"
fi

PYTHON="python3"
if [[ "$SKIP_VENV" -eq 0 ]]; then
    VENV_DIR="$REPO_ROOT/.venv-build"
    if [[ ! -d "$VENV_DIR" ]]; then
        echo "Creating build virtual environment at $VENV_DIR ..."
        "$PYTHON" -m venv "$VENV_DIR"
    fi
    PYTHON="$VENV_DIR/bin/python"
fi

echo "Using Python interpreter: $PYTHON"
"$PYTHON" --version

echo "Upgrading pip and installing build dependencies ..."
"$PYTHON" -m pip install --upgrade pip
"$PYTHON" -m pip install -r "$SCRIPT_DIR/requirements-build.txt"

echo "Running PyInstaller ..."
"$PYTHON" -m PyInstaller --noconfirm --clean "$SCRIPT_DIR/melanopsin_gui.spec"

shopt -s nullglob
apps=("$REPO_ROOT"/dist/MelanopsinModel-v*.app)
if [[ ${#apps[@]} -eq 0 ]]; then
    echo "Build finished but no .app was found in dist/." >&2
    exit 1
fi
APP="${apps[0]}"
ZIP_NAME="$(basename "${APP%.app}").zip"

echo "Zipping $(basename "$APP") -> $ZIP_NAME ..."
(
    cd "$REPO_ROOT/dist"
    ditto -c -k --keepParent "$(basename "$APP")" "$ZIP_NAME"
)

echo ""
echo "Build complete: $APP"
echo "Release zip:    $REPO_ROOT/dist/$ZIP_NAME"
echo "Place the .app inside a cloned copy of the repository so it can find"
echo "the data/ and myutils/ folders."

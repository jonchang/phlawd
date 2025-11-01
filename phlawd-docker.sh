#!/bin/bash
#
# Generic Docker runner script for phlawd programs
# This script is called via symlinks to determine which program to run

# Check if Docker is running
if ! docker info >/dev/null 2>&1; then
    echo "Error: Couldn't connect to Docker. (Run \`docker info\` for more details)"
    exit 1
fi

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_WAS_REMOTE=false

# Determine which Docker image to use
# Try local image first, then fall back to GitHub Container Registry
if docker image inspect phlawd:latest >/dev/null 2>&1; then
    IMAGE="phlawd:latest"
elif docker image inspect phlawd >/dev/null 2>&1; then
    IMAGE="phlawd"
else
    # Try to infer the GitHub repository from the current directory
    # Default to jonchang/phlawd if we can't determine it
    REPO_NAME="${GITHUB_REPOSITORY:-jonchang/phlawd}"
    IMAGE="ghcr.io/${REPO_NAME}:latest"
    CONTAINER_WAS_REMOTE=true
fi

# Determine which program to run based on $0
SYMLINK_NAME=$(basename "$0")

# Map symlink names to container binary names
case "$SYMLINK_NAME" in
    phlawd)
        PROGRAM_NAME="PHLAWD"
        ;;
    *)
        PROGRAM_NAME="$SYMLINK_NAME"
        ;;
esac

MOUNT_DIR="$(pwd)"

echo "Running ${SYMLINK_NAME} from container: ${IMAGE}"
echo "Mounting directory: ${MOUNT_DIR}"

if [ "$CONTAINER_WAS_REMOTE" = true ]; then
    echo "Using remote container from GitHub Container Registry"
    echo "To build a local image, read BUILD_DOCKER.md for instructions."
fi

echo ""

docker run --rm \
    -v "${MOUNT_DIR}:/workspace" \
    -w /workspace \
    "${IMAGE}" \
    "${PROGRAM_NAME}" "$@"


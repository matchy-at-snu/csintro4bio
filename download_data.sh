#!/bin/sh -e

hasCommand() {
    command -v "$1" >/dev/null 2>&1
}

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            FILENAME=$(basename "${OUTPUT}")
            DIR=$(dirname "${OUTPUT}")
            aria2c --max-connection-per-server="$(nproc)" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl; then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget; then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
    fail "No download tool found in PATH. Please install aria2c, curl or wget."
fi

# get project root
PROJECT_ROOT=$(git rev-parse --show-toplevel)

mkdir -p "$PROJECT_ROOT/data/hg38"
cd "$PROJECT_ROOT/data"

downloadFile "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz" "hg38.chromFa.tar.gz"

# Extract the downloaded file to the hg38 directory
tar -xzf hg38.chromFa.tar.gz -C hg38
rm hg38.chromFa.tar.gz

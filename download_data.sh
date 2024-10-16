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

curr_dir=$(pwd)

# Human
mkdir -p "$PROJECT_ROOT/data/hg38"
# chicken
mkdir -p "$PROJECT_ROOT/data/galGal3"
# C. elegans
mkdir -p "$PROJECT_ROOT/data/ce10"


cd "$PROJECT_ROOT/data"

# Human
downloadFile "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz" "hg38.chromFa.tar.gz"

# Chicken
downloadFile "https://hgdownload.soe.ucsc.edu/goldenPath/galGal3/bigZips/chromFa.tar.gz" "galGal3.chromFa.tar.gz"

# C. elegans
downloadFile "https://hgdownload.soe.ucsc.edu/goldenPath/ce10/bigZips/chromFa.tar.gz" "ce10.chromFa.tar.gz"

# Drosophila Melanogaster
downloadFile "https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz" "dm3.chromFa.tar.gz"

# Extract the downloaded file to the hg38 directory
echo "Extracting Human genome..."
pv hg38.chromFa.tar.gz | tar -xz -C hg38
echo "Extracting Chicken genome..."
pv galGal3.chromFa.tar.gz | tar -xz -C galGal3
echo "Extracting C. elegans genome..."
pv ce10.chromFa.tar.gz | tar -xz -C ce10
echo "Extracting Drosophila Melanogaster genome..."
pv dm3.chromFa.tar.gz | tar -xz -C dm3

rm hg38.chromFa.tar.gz
rm galGal3.chromFa.tar.gz
rm ce10.chromFa.tar.gz
rm dm3.chromFa.tar.gz

cd "$curr_dir"

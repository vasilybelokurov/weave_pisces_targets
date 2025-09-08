#!/bin/bash -e

# Pisces-specific WEAVE catalog generation script
# Uses local environment and database connection

# Get weave_tool directory and activate local environment  
TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$TOOL_DIR/weave_pisces_env/bin/activate"

# Set database host (use system PGUSER/PGPASSWORD)
export WSDB_HOST='wsdb.ast.cam.ac.uk'

# Add weave_chervin module to Python path
export PYTHONPATH="$TOOL_DIR/weave_fits_writeout/py:$PYTHONPATH"

# Pisces survey parameters
TRIMESTER='2025B2'
OFILE=WS2025B2-028-PISCES.fits
IPATH=$1

# Check if input file provided
if [ -z "$IPATH" ]; then
    echo "Usage: $0 INPUT_FITS_CAT"
    echo "Example: $0 /path/to/Pisces_enriched_021023.fits"
    exit 1
fi

# Check if input file exists
if [ ! -f "$IPATH" ]; then
    echo "Error: Input file $IPATH does not exist"
    exit 1
fi

echo "Creating WEAVE catalog for Pisces survey..."
echo "Input: $IPATH"
echo "Output: $OFILE"

# Run weave catalog generation tool directly with Python
python "$TOOL_DIR/weave_fits_writeout/bin/weave_chervin_make_catalog" \
           --output $OFILE \
           --progtemps='11331' \
           --trimester=$TRIMESTER  \
           --targcat=WS2025B2-028-PISCES.fits \
           --targsurvey=WS2025B2-028-PISCES \
           --gaia_version=3 \
           "$IPATH"

echo "WEAVE catalog created: $OFILE"
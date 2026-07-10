
INPUTS=/home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/year_foraccess.txt
CHUNK_DIR=/home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/chunks

mkdir -p "$CHUNK_DIR"

sed 's/[[:space:]]//g' "$INPUTS" | split -l 100 -d -a 4 - "$CHUNK_DIR/chunk_"

# chmod +x "$CHUNK_DIR/chunk_"
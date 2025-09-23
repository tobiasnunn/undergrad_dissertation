#!/bin/bash

EMAIL="0rtpjri2@anonaddy.me"
WORK_DIR="$HOME/source/01_inputs/04_fastas/_running_1500/"
LAST_COUNT_FILE="$HOME/.last_faa_count"

cd "$WORK_DIR"

# Get current count
current_count=$(ls -1 *_renamed.faa 2>/dev/null | wc -l)

# Get last reported count (or 0 if first run)
if [[ -f "$LAST_COUNT_FILE" ]]; then
    last_count=$(cat "$LAST_COUNT_FILE")
else
    last_count=0
fi

# Only send email if count has changed
if [[ $current_count -ne $last_count ]]; then
    echo "Prodigal job progress: $current_count files completed at $(date)" | \
        mail -s "Job Progress: $current_count files" "$EMAIL"
    echo "$current_count" > "$LAST_COUNT_FILE"
fi

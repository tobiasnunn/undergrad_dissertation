cd 01_inputs/04_fastas/_working/

# Get 5 random files and combine them
ls *_renamed.faa | shuf -n 5 | tee selected_files.txt | xargs cat > random_5_combined.faa

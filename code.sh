source=/share/pool/CompGenomVert/EpiAus/workflow/config
dest=/home/ychrysostomakis/barebones_assembly_dir/genome_assembly_github/config
rsync -av  --exclude '*/' --exclude '*/*/'  --update $source/ $dest/
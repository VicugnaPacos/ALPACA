# Instructions for the analysis of repeat coverage.
#
# Align reference to self with MUMmer
nucmer --maxmatch --nosimplify -p reference reference.fasta reference.fasta
# Convert output to tabular sorted by chromosome+position
show-coords -l -c -T -H reference.delta | sort -k12,12 -k1,2n > reference.coords
# Optionally list the repeats that will be processed 
list-repeats.pl reference.coords dummy 95 59
#
# Align assembly to reference with MUMmer
nucmer --maxmatch --nosimplify -p assembly reference.fasta assembly.fasta
# Convert output to tabular sorted by scaffold+position
show-coords -l -c -T -H assembly.delta | sort -k12,12 -k1,2n > assembly.coords
#
# Process the two alignment sets in tandem.
# Output lines contain: repeat_counter, repeat_unit_length, repeat_unit_separation, classification
repeat_content.pl reference.coords assembly.coords 95 95 > repeat.content
# Get subtotals per classification.
cut -d ' ' -f 4 repeat.content | sort | uniq -c
#


This is supposed to be a tool to make it easier for folks to just
chuck a modification into a genome.

Should take in some engineering command and sequence, then edit
the FASTA file and GFF file, and regenerate FAI for you.

Just thinking about homologous recombination from a PCR product.

===

Names:

- yuRecombinator
- Idealized Recombinator of GFF And FASTA (IRGAF)

===

Inputs:

- original FASTA and GFF files
- FASTA-esque file that has name and commands, currently RECOMBINE
- take the outside 30 and search the genome, find the location
- swap those sequences, bark and die noisily if weird

===

Procedure 

- read in all three files
- parse editing commands
- search function to find the sequence to recombine, and orientation
- function to swap the sequence, with right orientation, and then
  update the gff file and move everything appropriately

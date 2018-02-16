
# *H*elp *E*dit *N*ew *R*eference of *Y*east genome

This is a tool supposed to make it easier for folks to just
chuck a modification into a genome.

The aim is to take a construct in genbank format, and recombine
that into the reference yeast genome by recombining based on left and
right arm perfect homology. It then generates a new indexed FASTA
file of the genome and a GFF annotation, with the genbank annotations
carried over as approximate GFF3 annotations.

## ToDo

- Make it take options, for a more reasonable interface.
- Make it use the reverse (swap the indicies). Currently, it just
  searches the plus strand.
- Read in Genbank format for the construct !!!
- Figure out a better way to call the `faidx` program (or generate
  it in script???)

This is supposed to be a tool to make it easier for folks to just
chuck a modification into a genome.

Basically, it's gonna take two GFF3 files, a genome and a construct.
It looks for perfect matches to start extending from, figures all
the perfect homology, then just sticks the new construct sequence
in between the two homology bits on the genome.

===

TODO

- make it search all chromosomes (time thing, not for debugging)
- make it use the reverse (swap the indicies)
- make options and crap
- figure better way to call faidx
- call IGV?

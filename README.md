
# **H**elp **U**pdate **GFF** **F**iles

**HUGFFF**

This is a tool supposed to make it easier for folks to just
chuck a modification into a genome.

You need to provide a GFF 3 with sequence in it, and then a genbank
formatted construct. Then it takes homology off each end, finds 
seeds, extends those, then splices together the genome based on those
coordinates.
Nothing fancy. It's pretty much just BioPython munging.
Thanks BioPython!

There's a GFF3 readily available for download from SGD, so it's easy
for yeast people.

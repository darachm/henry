
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

# Usage-ish

As seen when you run without arguments.
The only required ones are `--referenceGFF` which points toward the reference
GFF3 file, `--constructGenBank` is the genbank file of the thing you are 
transforming in, and `--output_base` is the basename of the new files.

    usage: recombinator.py [-h] --referenceGFF REFERENCEGFF --constructGenBank
                           CONSTRUCTGENBANK [--homology_start HOMOLOGY_START]
                           [--homology_end HOMOLOGY_END]
                           [--homology_start_5 HOMOLOGY_START_5]
                           [--homology_end_5 HOMOLOGY_END_5]
                           [--homology_start_3 HOMOLOGY_START_3]
                           [--homology_end_3 HOMOLOGY_END_3]
                           [--mismatches MISMATCHES] [--mismatches_5 MISMATCHES_5]
                           [--mismatches_3 MISMATCHES_3] --output_base OUTPUT_BASE

Other parameters deal with tolerable mismatches and adjustable homology on each
end of the construct you're trying to transform in.

It's not going to find homology in the middle of the construct, just the very
ends of it.

# Dependencies

You're going to need `BioPython` and `BCBio`. So this can be done with `pip`,
so `pip3 install BioPython` and `pip3 install bcbio-gff`.

If you don't have `pip`, ... ?

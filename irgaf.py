#!/usr/bin/python3
# This thing is going to take a PCR product, transform it into a
# genome using 30bp homology, and then update the fasta and GFF to 
# reflect that.

from Bio import SeqIO
from Bio import motifs
from Bio import Alphabet
from BCBio import GFF

# Configuration, define these vars
#genome_fasta = "S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa"
gff_file = "S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff"
construct_file = "construct.fa"
# Later, should be setup as arguments

#the_genome = list(SeqIO.parse(genome_fasta,"fasta"))
the_gff = list(GFF.parse(gff_file))
construct = list(SeqIO.parse(construct_file,"fasta",Alphabet.DNAAlphabet()))

#for record in the_gff:
#  print(record.features)

def find_recombination_site(gff,construct,search_length=24):
  construct_end_left= motifs.create([construct[:search_length].seq],alphabet=construct.seq.alphabet)
  construct_end_right = motifs.create([construct[-search_length:].seq],alphabet=construct.seq.alphabet)
#
#  for chromosome in gff:
  chromosome = gff[10]
  chromosome.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
#
#  search_left = construct_end_left.counts.normalize(pseudocounts=0.1).log_odds().search(chromosome.seq,threshold=30)
  left_match = list(
    construct_end_left.instances.search(chromosome.seq)
    ,construct_end_left.instances.search(chromosome.seq)
    )[0]
  right_match = list(construct_end_right.instances.search(chromosome.seq))[0]
  return(left_match+right_match)
#  return(construct.seq+"\n"+construct_end_left.seq+"\n"+construct_end_right.seq)
  


print(find_recombination_site(the_gff,construct[0]))



#for position, score in pssm.search(test_seq, threshold=3.0):
#     print("Position %d: score = %5.3f" % (position, score))
#pos is located at test_seq[pos:pos+len(m)] 


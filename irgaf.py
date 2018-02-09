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

the_gff = list(GFF.parse(gff_file))
construct = list(SeqIO.parse(construct_file,"fasta"))

#
#
#
#
#
#
# rewrite to just shrink in the homology to the first base of change
#
#
#
#
#
#

def find_recombination_sites(gff,construct,search_length=24):
  construct_end_left= motifs.create([construct[:search_length].seq])
  construct_end_right = motifs.create([construct[-search_length:].seq])
  for chromosome_number,chromosome in enumerate(gff):
    if chromosome_number != 10:
      continue
    chromosome.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
#  search_left = construct_end_left.counts.normalize(pseudocounts=0.1).log_odds().search(chromosome.seq,threshold=30)
    matches = []
    for left_seed in construct_end_left.instances.search(chromosome.seq):
      for right_seed in construct_end_right.instances.search(
          chromosome.seq[left_seed[0]:]
          ):
        a_match = (chromosome_number
          ,left_seed[0]
          ,left_seed[0]+right_seed[0]+search_length
          ,"+"
          )
        matches.append(a_match)
    rev_chromosome =  chromosome.reverse_complement()
    for left_seed in construct_end_left.instances.search(rev_chromosome.seq):
      for right_seed in construct_end_right.instances.search(
          rev_chromosome.seq[left_seed[0]:]
          ):
        a_match = (chromosome_number
          ,len(rev_chromosome.seq)-left_seed[0]
          ,len(rev_chromosome.seq)-
            (left_seed[0]+right_seed[0]+search_length)
          ,"-"
          )
        matches.append(a_match)
  length_construct = len(construct)
  print("Construct is "+str(length_construct)+"bp long.")
  print("I found homologies at spans:")
  for each_match in matches:
    print("\t"+str(each_match[2]-each_match[1])+"bp long on "+
      "chromosome "+gff[each_match[0]].id+", from "+
      str(each_match[1])+" to "+str(each_match[2]))
  return(matches)

def recombine_at_match(gff,construct,this_match):
  chromosome_number = this_match[0]
  match_start = this_match[1]
  match_end = this_match[2]
  match_length = match_end-match_start
  match_strand = this_match[3]
  construct_length = len(construct)
  shift_length = construct_length-match_length


  
some_matches = find_recombination_sites(the_gff,construct[0])
#print(recombine_at_match(the_gff,construct[0],(10,513905,513985,"+")))


#!/usr/bin/python3
# This thing is going to take a PCR product, transform it into a
# genome using 30bp homology, and then update the fasta and GFF to 
# reflect that.

from Bio import SeqIO
from Bio import motifs
from Bio import Alphabet
from BCBio import GFF
from subprocess import call

# Configuration, define these vars
# Later, should be setup as arguments
gff_file = "S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff"
construct_file = "BarcodeLibrary1.gff"
homology_region = [0,24] # how bit of a slice to take off the ends?


output_base = "test"



the_gff   = list(GFF.parse(gff_file))
#construct = list(SeqIO.parse(construct_file,"fasta"))[0]
construct = list(GFF.parse(construct_file))[0]





def find_recombination_sites(gff,construct
  ,left_homology_slice,right_homology_slice
  ):
  
  # For each chromosome, without a dict lookup:
  for chromosome_number,chromosome in enumerate(gff):
    
    # For debugging, for speed
#    if chromosome_number != 1:
#      continue
    
    # Because alphabets are hard to inherit I suppose
    chromosome.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
    
    # List to hold all matches
    matches = []
    
    # We start looking from the plus strand, for each match on the
    # left
    for each_result in (
      list_of_slices_for_recombine(chromosome
        ,construct,left_homology_slice,right_homology_slice
        )
      ):
      matches.append([chromosome_number,"+",each_result])
    
    # Then we flip around the rev_chromosome and look from the other end
    rev_chromosome =  chromosome.reverse_complement()
    for each_result in (
      list_of_slices_for_recombine(rev_chromosome
        ,construct,left_homology_slice,right_homology_slice
        )
      ):
# NEED TO FLIP AROUND COORDS TO FORWARD!!!
      matches.append([chromosome_number,"-",each_result])
  
  return(matches)






def list_of_slices_for_recombine(reference,construct
    ,left_slice,right_slice
    ):
  
  recombines = []
  
  for left_match in (
    extended_matches(reference,construct,left_slice)
    ):
    
    for right_match in (
      extended_matches(reference,construct,right_slice)
      ):
      
      recombines.append([left_match,right_match])
  
  return(recombines)





def extended_matches(chromosome,construct,seed_slice):
  
  # We start in from the left of the sequence, and try to find
  # perfect match, extend it a bit, then return start and stop
  
  seed_motif = motifs.create([construct.seq[slice(*seed_slice)]])
  
  matches = []
  
    # This is for fuzzy matching 
    #  search_left = (construct_end_left.counts
    #    .normalize(pseudocounts=0.1).log_odds()
    #    .search(chromosome.seq,threshold=30)
  
  for seed_pos,seed_seq in (
    seed_motif.instances.search(chromosome.seq)
    ):
  
    match_start = seed_pos
    match_end   = seed_pos+len(seed_seq)
  
    try:
      while (construct.seq[match_start-seed_pos+seed_slice[0]] == 
        chromosome.seq[match_start]
        ):
        match_start -= 1
      else:
        match_start += 1
    except:
      match_start += 1
  
    try:
      while (construct.seq[match_end-seed_pos+seed_slice[0]] == 
        chromosome.seq[match_end]
        ):
        match_end += 1
      else:
        match_end -= 1
    except:
      match_end -= 1
    
    # These are where things match
    reference_match = [match_start,match_end]
    construct_match = [match_start-seed_pos+seed_slice[0]
      ,match_end-seed_pos+seed_slice[0]
      ]
    
    # Here, I've got to flip it around because indicies are weird
    # in python with regards to the end...
#DELETE?
#    for i,j in enumerate(construct_match):
#      if construct_match[i] < 0:
#        construct_match[i] = len(construct.seq)-construct_match[i]+1
    
    matches.append([reference_match,construct_match])
  
  return(matches)





def recombine_at_match(gff,construct,this_match):
  
  chromosome = gff[this_match[0]]
  
  excise_start = this_match[2][0][0][1]
  excise_end = this_match[2][1][0][0]
  excise_length = excise_end - excise_start
  
  integrate_start = this_match[2][0][1][1]
  integrate_end = this_match[2][1][1][0]
  integrate_length = integrate_end - integrate_start
  
  shift_length = integrate_length-excise_length
  
  gff[this_match[0]] = (
    chromosome[:excise_start] +
    construct[integrate_start:integrate_end] +
    chromosome[excise_end:]
    )
  gff[this_match[0]].id = chromosome.id
  
  return(gff)


#" part of ~/.vimrc
#" highlight tabs and trailing spaces
#set listchars=tab:>-,trail:-
#set list  


if __name__ == "__main__":
  
  some_matches = find_recombination_sites(the_gff
    ,construct
    ,[3,24],[len(construct.seq)-24+1,len(construct.seq)-3+1]
    )

  print("I may have found a match ? ")
  print(some_matches)

  for each_match in some_matches:
    the_gff = recombine_at_match(the_gff,construct,each_match)

  with open(output_base+".gff", "w") as out_gff, open(output_base+".fa", "w") as out_fasta:
    GFF.write(the_gff,out_gff)
    out_fasta.write("")

  with open(output_base+".gff", "a") as out_gff,open(output_base+".fa", "a") as out_fasta:
    out_gff.write("##FASTA\n")
    for i,chromosome in enumerate(the_gff):
      out_gff.write(">"+chromosome.id+"\n")
      out_gff.write(str(chromosome.seq)+"\n")
      out_fasta.write(">"+chromosome.id+"\n")
      out_fasta.write(str(chromosome.seq)+"\n")
  
  call("samtools faidx "+output_base+".fa",shell=True)


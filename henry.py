#!/usr/bin/python3
# This thing is going to take a PCR product, transform it into a
# genome using 30bp homology, and then update the fasta and GFF to 
# reflect that.

import argparse
import copy
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqFeature
from Bio import SeqIO
from Bio import motifs
from Bio import Alphabet
from BCBio import GFF
from subprocess import call


def find_recombination_sites( reference, construct, 
        left_homology_slice, right_homology_slice ):
    
    # List to hold all matches
    matches = []

    # For each chromosome, without a dict lookup:
    for chromosome_number, chromosome in enumerate(reference):
    
        print("Looking on chromosome "+chromosome.id)
        
        # For debugging, for speed
        if chromosome_number != 1:
            continue

        # We start looking from the plus strand, for each match on the
        # left
        for each_result in (
            list_of_slices_for_recombine(chromosome, construct,
                left_homology_slice,right_homology_slice
                )
            ):
            matches.append([chromosome_number,"+",each_result])
            print("\tI may have found a match ? ")
            print("\t+ "+str(each_result))
        
        # Then we flip around the rev_chromosome and look from the other end
        rev_chromosome = chromosome.reverse_complement()
        for each_result in (
            list_of_slices_for_recombine(rev_chromosome, construct,
                left_homology_slice,right_homology_slice
                )
            ):
            matches.append([chromosome_number,"-",each_result])
# NEED TO FLIP AROUND COORDS TO FORWARD!!!
            print("\tI may have found a match ? ")
            print("\t- "+str(each_result))
    
    return(matches)


def list_of_slices_for_recombine(reference, construct
        ,left_slice,right_slice
        ):
    
    recombines = []
    
    for left_match in (
        extended_matches(reference,construct,left_slice,args.mismatches)
        ):
        
        for right_match in (
            extended_matches(reference,construct,right_slice,args.mismatches)
            ):
            
            recombines.append([left_match,right_match])
    
    return(recombines)

def search_for_motif(target,query,mode="perfect",mismatches=0):

    query = Seq.Seq(str(query.seq),
        alphabet=Alphabet.IUPAC.IUPACUnambiguousDNA())
    seed_motif = motifs.create([query])

    target = Seq.Seq(str(target.seq),
        alphabet=Alphabet.IUPAC.IUPACUnambiguousDNA())

    if mode is "fuzzy":
    
        pssm = ( seed_motif.counts
                .normalize(pseudocounts=0.1).log_odds() )
        
        pssm_length = len(query)
        pssm_max = pssm.max/pssm_length
        pssm_min = pssm.min/pssm_length
    
        searches = pssm.search(target,
            threshold=-0.001+(pssm_length-int(mismatches))*pssm_max+
                int(mismatches)*pssm_min
            )

    elif mode is "perfect":

        searches = seed_motif.instances.search(target)

    return(searches)


def extended_matches(chromosome,construct,seed_slice,mismatches=0):
    
    # We start in from the left of the sequence, and try to find
    # perfect match, extend it a bit, then return start and stop

    query_seq = construct[slice(*seed_slice)]
    
    searches = search_for_motif(chromosome,query_seq,mode="perfect")

    matches = []

    for seed_pos,seed in searches :
    
        match_start = seed_pos
        match_end   = seed_pos+len(seed)

        try:
            while (construct.seq[match_start-seed_pos+seed_slice[0]] == 
                    chromosome.seq[match_start]
                    ):
                if (match_start-seed_pos+seed_slice[0]) <= 0:
                    break
                match_start -= 1
            else:
                match_start += 1
        except:
            match_start += 1
    
        try:
            while (construct.seq[match_end-seed_pos+seed_slice[0]] == 
                chromosome.seq[match_end]
                ):
                if (match_end-seed_pos+seed_slice[0]) > len(construct.seq):
                    break
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
        
        matches.append([reference_match,construct_match])
    
    return(matches)


def recombine_at_match(gff,construct,this_match):

    chromosome = gff[this_match[0]]
    
    if this_match[1] == "+":
        excise_start = this_match[2][0][0][1]
        excise_end = this_match[2][1][0][0]
    else:
        excise_end = len(chromosome) - this_match[2][0][0][1]
        excise_start = len(chromosome) - this_match[2][1][0][0]
    excise_length = excise_end - excise_start
    
    if this_match[1] == "+":
        integrate_start = this_match[2][0][1][1]
        integrate_end = this_match[2][1][1][0]
    else:
        integrate_end = len(construct) - this_match[2][0][1][1]
        integrate_start = len(construct) - this_match[2][1][1][0]
    integrate_length = integrate_end - integrate_start

    print("Recombining reference, splicing on chromosome #"+
        str(this_match[0]+1)+" after "+str(excise_start)+
        ", starting at position "+str(integrate_start)+
        " in the construct through to "+str(integrate_end)+
        ", then back to "+str(excise_end)+" in the reference")
    
    gff[this_match[0]] = (
        chromosome[:excise_start] +
        construct[integrate_start:integrate_end] +
        chromosome[excise_end:]
        )
    gff[this_match[0]].id = chromosome.id
    
    return(gff)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=""+
        "slapChopper.py")
    parser.add_argument("--referenceGFF",
        default="S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff")
    parser.add_argument("--constructGenBank",
        default="BarcodeLibrary1.gb")
    parser.add_argument("--homology_start",   default=0)
    parser.add_argument("--homology_end",     default=20)
    parser.add_argument("--homology_start_5", default=None)
    parser.add_argument("--homology_end_5",   default=None)
    parser.add_argument("--homology_start_3", default=None)
    parser.add_argument("--homology_end_3",   default=None)
    parser.add_argument("--mismatches",       default=0)
    parser.add_argument("--mismatches_5",     default=None)
    parser.add_argument("--mismatches_3",     default=None)
    parser.add_argument("--output_base",      default="output")
    args=parser.parse_args()

    if args.homology_start_5 is None:
        args.homology_start_5 = args.homology_start
    if args.homology_end_5 is None:
        args.homology_end_5 = args.homology_end
    if args.homology_start_3 is None:
        args.homology_start_3 = args.homology_start
    if args.homology_end_3 is None:
        args.homology_end_3 = args.homology_end

    if args.mismatches_5 is None:
        args.mismatches_5 = args.mismatches
    if args.mismatches_3 is None:
        args.mismatches_3 = args.mismatches

    print("===")
    for key, value in vars(args).items():
        print("Argument\t"+key+"\tis\t"+str(value))
    print("===")

    the_reference = list(GFF.parse(args.referenceGFF))
    for chromosome in the_reference:
        # Because alphabets are hard to inherit I suppose
        chromosome.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()

    constructs = list(SeqIO.parse(args.constructGenBank, "genbank"))

    for each_construct in constructs:

####TESTING
        each_construct = each_construct.reverse_complement()
####TESTING
    
        features_copy = list(each_construct.features)

        for i, each_feature in enumerate(features_copy):

            qualifiers = dict(each_feature.qualifiers)
            for each_key in list(qualifiers.keys()):
                qualifiers[each_key.replace("ApEinfo_","")] = \
                    qualifiers[each_key]

            if "color" not in qualifiers.keys():

                try:
                    qualifiers["color"] = qualifiers["fwdcolor"]
                except:
                    pass

            each_construct.features[i] = SeqFeature.SeqFeature(
                location=each_construct.features[i].location,
                type=each_construct.features[i].type,
                strand=each_construct.features[i].strand,
                ref=each_construct.features[i].ref,
                ref_db=each_construct.features[i].ref_db,
                qualifiers=qualifiers
                )

        each_construct.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
    
        some_matches = find_recombination_sites(
            the_reference, each_construct,
            [args.homology_start_5, args.homology_end_5],
            [len(each_construct.seq)-args.homology_end_3+1, 
                len(each_construct.seq)-args.homology_start_3+1] 
            )

        print()
        print("I found "+str(len(some_matches))+" matches.")
        print()
    
        for each_match in some_matches:
            the_reference = recombine_at_match(the_reference,
                    each_construct,each_match)

    with open(args.output_base+".gff", "w") as out_gff, \
            open(args.output_base+".fa", "w") as out_fasta:
        GFF.write(the_reference,out_gff)
        out_fasta.write("")

    with open(args.output_base+".gff", "a") as out_gff, \
            open(args.output_base+".fa", "a") as out_fasta:
        out_gff.write("##FASTA\n")
        for i,chromosome in enumerate(the_reference):
            out_gff.write(">"+chromosome.id+"\n")
            out_gff.write(str(chromosome.seq)+"\n")
            out_fasta.write(">"+chromosome.id+"\n")
            out_fasta.write(str(chromosome.seq)+"\n")
    
    try:
        call("samtools faidx "+args.output_base+".fa",shell=True)
    except:
        raise("well you don't have samtools configured for use")



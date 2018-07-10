#!/usr/bin/python3

# This script takes in a reference GFF file, a genbank formatted
# construct, then searches around for homology to recombine it in
# on one chromosome, applies those changes and returns a GFF file.

import argparse
import copy
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqFeature
from Bio import SeqIO
from Bio import motifs
from Bio import Alphabet
from BCBio import GFF # This BCBio is installable as bcbio-gff in pip
from subprocess import call

def find_recombination_sites( reference, construct, 
        left_homology_slice, right_homology_slice, 
        mismatches_5=0, mismatches_3=0 ):

    # This looks for recombination sites, given a reference and
    # construct sequence, and also slices of how much to use and
    # what mismatches are tolerated (to trigger fuzzy matching)
    
    # List to hold all matches
    matches = []

    # For each chromosome
    for chromosome_number, chromosome in enumerate(reference):
    
        # Note that these are biopython Seqs or SeqRecords now
        print("Looking on chromosome "+chromosome.id)

        # For debugging you can uncomment this to speed things up...
#        if chromosome_number != 1:
#            continue
        
        # We start looking from the plus strand, for each match on the
        # left
        for each_result in (
            list_of_slices_for_recombine(chromosome, construct,
                left_homology_slice,right_homology_slice, 
                mismatches_5=mismatches_5, mismatches_3=mismatches_3
                )
            ):
            # If we get one, we keep it
            matches.append([chromosome_number,"+",each_result])
            # and report it
            print("\tI may have found a match ? ")
            print("\t+ "+str(each_result))
        
        # Then we flip around the rev_chromosome and look from the other end
        rev_chromosome = chromosome.reverse_complement()
        for each_result in (
            list_of_slices_for_recombine(rev_chromosome, construct,
                left_homology_slice,right_homology_slice, 
                mismatches_5=mismatches_5, mismatches_3=mismatches_3
                )
            ):
            matches.append([chromosome_number,"-",each_result])
            print("\tI may have found a match ? ")
            print("\t- "+str(each_result))
    
    # Then we return now all the matches from left and right on that
    # reference genome
    return(matches)


def list_of_slices_for_recombine(reference, construct
        ,left_slice, right_slice,
        mismatches_5=0, mismatches_3=0,
        ):

    # Here we want to return the slices for recombination, just for
    # a reference and construct, and a certain slice off the ends
    
    recombines = []
    
    # We start by looking from the left, and when we find an 
    # extended match
    for left_match in (
        extended_matches(reference, 0,
            construct, left_slice, mismatches_5 )
        ):

        # Then we continue from that match along to the right,
        # looking for a right match to finish recombining at.
        for right_match in (
            extended_matches(reference, left_match[0][0],
                construct, right_slice, mismatches_3 )
            ):

            # Note there's some inefficiencies here. If I match twice
            # coming from the left, and twice looking to the right,
            # then I have to re-match the two right ones for each
            # left match, if that makes sense. Would be more 
            # efficient to find all matches, index them, then feed
            # them in as needed. But, this was easier to write.
            
            # If we find a match, append
            recombines.append([left_match,right_match])
    
    # Return all recombination slices
    return(recombines)


def extended_matches(chromosome,offset,construct,seed_slice,mismatches=0):
    
    # We start in from the left of the sequence, and try to find
    # perfect match, extend it a bit, then return start and stop

    # We've got a query seed
    query_seq = construct[slice(*seed_slice)]

    # And the chromosome can be partial, especially if we're looking
    # for that second match for recombination
    chromosome = chromosome[offset:]
    
    # Look for the motif
    searches = search_for_motif(chromosome, query_seq, mismatches=mismatches)

    matches = []

    for seed_pos, seed in searches :
    
        # For each seed that sticks, we get the positions
        match_start = seed_pos
        match_end   = seed_pos+len(seed)

        # Then we try to walk along and test if we can resect the
        # start of it back on both construct and genome, but only
        # where there's perfect matching. If there isn't and it's
        # not biologically relevant, than it'll just recombine here.
        # If you've got a SNP way outside your altered construct
        # that's annotated, then that just means that you'll lose the
        # native annotations in the genome.
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
    
        # We want to resect it forwards so that we can maintain
        # annotations that are mainly in the reference GFF, so the
        # resect on the forward doesn't really matter at all.
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
        
        # These are where things match, so let's build it nicely
        reference_match = [match_start+offset,match_end+offset]
        construct_match = [match_start-seed_pos+seed_slice[0]
            ,match_end-seed_pos+seed_slice[0]
            ]
        
        # Stick it on our list of matches
        matches.append([reference_match,construct_match])
    
    # And return them all
    return(matches)


def search_for_motif(target,query,mismatches=0):

    # Here we're just searching for a motif 

    # We copy the query sequence we're searching for so that we can
    # make sure the alphabet is Unambiguous
    query = Seq.Seq(str(query.seq),
        alphabet=Alphabet.IUPAC.IUPACUnambiguousDNA())

    # Then we make a motif
    seed_motif = motifs.create([query])

    # We do the same Unambiguous alphabet on the target, 
    # the reference genome
    target = Seq.Seq(str(target.seq),
        alphabet=Alphabet.IUPAC.IUPACUnambiguousDNA())

    if mismatches == 0:

        # This is perfect matches only, super faster
        searches = seed_motif.instances.search(target)

    elif mismatches > 0:
    
        # Otherwise, we make a position scoring matrix, and convert
        # to log odds so we can scan it along the genome
        pssm = ( seed_motif.counts
                .normalize(pseudocounts=0.1).log_odds() )
        
        # Some little helper vars
        pssm_length = len(query)
        pssm_max = pssm.max/pssm_length
        pssm_min = pssm.min/pssm_length
    
        # Then we search with a threshold of mismatches tolerated
        searches = pssm.search(target,
            threshold=-0.001+(pssm_length-int(mismatches))*pssm_max+
                int(mismatches)*pssm_min
            )

        # This is sooooo slow. A real tool should have some indexing
        # to quicken the search/alignment. Maybe some sort of bash
        # script using a BWA alignment to seed, then extend and
        # splice using a BioPython script like this.

    return(searches)


def recombine_at_match(gff,construct,this_match):

    # Here we regenerate the gff as we need it

    # Here, we pick out the chromosome we had a hit on
    chromosome = gff[this_match[0]]
    
    # Then we identify the start and stop coordinates of the match
    # along the genome. Sorry about the depths of indexing, it's
    # not transparent and you'll just have to read the function
    # that returned the matches.
    if this_match[1] == "+":
        excise_start = this_match[2][0][0][1]
        excise_end = this_match[2][1][0][0]
    else:
        excise_end = len(chromosome) - this_match[2][0][0][1]
        excise_start = len(chromosome) - this_match[2][1][0][0]
    
    # Then we identify the start and stop coordinates of the match
    # along the construct
    integrate_start = this_match[2][0][1][1]
    integrate_end = this_match[2][1][1][0]
    if this_match[1] == "+":
        the_chunk_going_in = construct[integrate_start:integrate_end]
    else:
        the_chunk_going_in = \
            construct[integrate_start:integrate_end]\
            .reverse_complement()

    # We report what's up
    print("Recombining reference, splicing on chromosome #"+
        str(this_match[0]+1)+" after "+str(excise_start)+
        ", starting at position "+str(integrate_start)+
        " in the construct through to "+str(integrate_end)+
        ", then back to "+str(excise_end)+" in the reference.")
    
    # Then we just splice the reference chromosome from the end of
    # the homology, then the construct from the end of the left
    # homology to the start of the right homology, then from the 
    # start of the homology in the genome to the end of the 
    # chromosome.
    gff[this_match[0]] = (
        chromosome[:excise_start] +
        the_chunk_going_in +
        chromosome[excise_end:]
        )
    gff[this_match[0]].id = chromosome.id

    return(gff)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="recombinator.py"+
        "Help Update GFF Files (HUGFFF)."+
        "It's designed to take a construct in as "+
        "a GenBank formatted sequence record (so ApE does that as "+
        "an export), and a reference GFF (you can get it from SGD "+
        "downloads), search for homology, and splice together a "+
        "reference. YMMV, check the work for little indels to make"+
        "sure the indexing is right. And re-implement if you want"+
        "the fuzzy matching to run in this lifetime."+
        ""+
        "You can specify the starts and stops for each end of the"+
        "construct with the arguments.")
    parser.add_argument("--referenceGFF",     required=True)
    parser.add_argument("--constructGenBank", required=True)
    parser.add_argument("--homology_start",   default=0)
    parser.add_argument("--homology_end",     default=20)
    parser.add_argument("--homology_start_5", default=None)
    parser.add_argument("--homology_end_5",   default=None)
    parser.add_argument("--homology_start_3", default=None)
    parser.add_argument("--homology_end_3",   default=None)
    parser.add_argument("--mismatches",       default=0)
    parser.add_argument("--mismatches_5",     default=None)
    parser.add_argument("--mismatches_3",     default=None)
    parser.add_argument("--output_base",      required=True)
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
    print("SETTINGS:")
    for key, value in vars(args).items():
        print("Argument : "+key.ljust(16)+"\tis\t"+str(value))

    if int(args.mismatches) > 0:
        print()
        print("FYI, you set mismatches > 0. That is slooooooooow. "+
            "If you can get by with using 0 mismatches, please do "+
            "so.") 
        print()

    # We want a GFF 3 with the sequence information at the bottom,
    # because that's our source of the chromosome info.
    the_reference = list(GFF.parse(args.referenceGFF))
    for chromosome in the_reference:
        # Because alphabets are hard to inherit I suppose
        chromosome.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
    # You could modify this to take it from a different source, but
    # then you'd want to modify this to make a BioPython SeqRecord
    # with all the annotation in it. Your choice. The GFF3 is readily
    # downloadable from SGD, and I really like the atomic format of
    # it, I think the two together is real smart.

    # Here we read in a GenBank construct
    try:
        constructs = list(SeqIO.parse(args.constructGenBank, "genbank"))
    except:
        raise("Well, I couldn't read that construct in. "+
            "Maybe that's not real genbank format?")

    # Here we clean up that construct or constructs (you can put
    # multiple in a GenBank maybe?). This is primarily for IGV to
    # make the colors pretty, by renaming them as color. Mainly for
    # ApE outputs.
    for each_construct in constructs:

        # Make a copy for iterating
        features_copy = list(each_construct.features)

        for i, each_feature in enumerate(features_copy):

            # Make a copy for iterating
            qualifiers = dict(each_feature.qualifiers)

            for each_key in list(qualifiers.keys()):

                # We're just cutting out the ApEinfo_ to make it
                # more of a general genbank format
                qualifiers[each_key.replace("ApEinfo_","")] = \
                    qualifiers[each_key]

            # If it's not defined, then set it as fwdcolor
            if "color" not in qualifiers.keys():

                try:
                    qualifiers["color"] = qualifiers["fwdcolor"]
                except:
                    pass

            # Redefine the feature with the new qualifiers
            each_construct.features[i] = SeqFeature.SeqFeature(
                location=each_construct.features[i].location,
                type=each_construct.features[i].type,
                strand=each_construct.features[i].strand,
                ref=each_construct.features[i].ref,
                ref_db=each_construct.features[i].ref_db,
                qualifiers=qualifiers
                )

        each_construct.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
    
        # Go find some matches
        some_matches = find_recombination_sites(
            the_reference, each_construct,
            [int(args.homology_start_5), int(args.homology_end_5)],
            [len(each_construct.seq)-int(args.homology_end_3)+1, 
                len(each_construct.seq)-int(args.homology_start_3)+1],
            mismatches_5=int(args.mismatches_5),
            mismatches_3=int(args.mismatches_3)
            )

        print()
        print("I found "+str(len(some_matches))+" matches.")
        print()
    
        # For each match in series, edit the genome. Note that
        # each construct gets edited after each other, so if you
        # can put multiple in a GenBank (can you???) then they could
        # edit each other's edits.
        for each_match in some_matches:
            the_reference = recombine_at_match(the_reference,
                    each_construct,each_match)

    # Finally, write out the recombined files, both as a GFF and
    # as a fasta, and index it for you so you can IGV easily.
    with open(args.output_base+".gff", "w") as out_gff, \
            open(args.output_base+".fa", "w") as out_fasta:
        GFF.write(the_reference,out_gff)
        out_fasta.write("")
    with open(args.output_base+".gff", "a") as out_gff, \
            open(args.output_base+".fa", "a") as out_fasta, \
            open(args.output_base+".fa.fai", "w") as out_fasta_index:
        out_gff.write("##FASTA\n")
        index_offset = 0
        for i,chromosome in enumerate(the_reference):
            out_gff.write(">"+chromosome.id+"\n")
            out_gff.write(str(chromosome.seq)+"\n")
            out_fasta.write(">"+chromosome.id+"\n")
            out_fasta.write(str(chromosome.seq)+"\n")
            out_fasta_index.write(chromosome.id+"\t"+
                str(2+len(chromosome.seq))+"\t"+
                str(2+len(chromosome.id)+index_offset)+"\t"+
                str(2+len(chromosome.id)+len(chromosome.seq))+"\t"+
                str(2+len(chromosome.id)+len(chromosome.seq)+1)+"\n"
                )
            index_offset += 3 + len(chromosome.id) + len(chromosome.seq)


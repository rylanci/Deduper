#!/usr/bin/env python3

"""Deduper Code by Ryan Lancione"""

#Packages
import argparse
import re
import os


def get_args():
    parser = argparse.ArgumentParser(description="A program to remove PCR duplicates from sam files")
    parser.add_argument("-i", help="Input sam file", required = True)
    parser.add_argument("-j", help="Output sam file", required = True)
    parser.add_argument("-u", help="file with unique molecular identifiers or set to 'random'")
#    parser.add_argument("-s", help="If input sam sorting desired, state 'sort'", required = False)
    parser.add_argument("-p", help = "If paired end sequencing data state 'True'", required = False)
    return parser.parse_args()

args = get_args()

#If paired end mode selected print an error and exit the program.
if args.p == "True":
    print("Paired end capabilities not yet implemented")
    exit()



#Input sam file
input = open(args.i, "r")

'''
if args.s == "sort":
    sam_to_bam = f"samtools view -S -b {input} ">" temp.bam"
    sort_bam = "samtools sort temp.bam - O sam -o sorted.sam"
    os.system(sam_to_bam)
    os.system(sort_bam)
'''
#Output sam file
output = open(args.j, "w")

#Umi_list. implement in argparse
umi_file = open(args.u, 'r')
umi_list = umi_file.read().splitlines()


#Funtions
def mapped(bitflag):
    '''
    Check 0x4 of the bitflag to tell if the read was mapped
    Rewrite to take bitflag as param, but make save line split for final loop
    '''
    if ((bitflag & 4) != 4):
	#mapped_true
        return True
    else:
	#mapped_false
        return False


def fwORrv(bitflag):
    '''
    Reads bit flag to ask whether the read aligned to the forward or reverse strand.
    '''
    if ((bitflag & 16) != 16):
	    #Forward seq aligning
        return True
    else:
	   #Reverse seq aligning
       return False


#Our set that will contain all the unique records
#It will reset for each new chromosome to manage memory
uniq_recs = set()
previous_chrom = 0


#***Main Loop***
for line in input:

    #Write header to file
    if line[0] == "@":
        output.write(line)

    #Start Deduping!
    elif line[0] != "@":
        #Store our important record attributes
        qname = line.split()[0]
        bitflag = int(line.split()[1])
        lpos = line.split()[3]
        CIGAR = line.split()[5]
        chrom = qname.split(":")[3]
        umi = qname.split(":")[7]

        #Flush our set if the chromosome changes
        if chrom != previous_chrom:
            uniq_recs.clear()
        previous_chrom = chrom


        #Check if umi in umi set and mapped
        # if so parse bitflag to tell whether read is fw or rv
        #Adjust left align coords approprately
        if umi in umi_list or args.u == "random" and mapped(bitflag) == True:

            #If forward aligning subtract for int softclipped
            cigar_tups = re.findall(r'([0-9]+)([A-Z]{1})', CIGAR)
            if fwORrv(bitflag):
                if "S" in cigar_tups[0][1]:
                    adj_lpos = int(lpos) - int(cigar_tups[0][0])
                else:
                    adj_lpos = int(lpos)

            #If reverse calculate the 5' position
            if fwORrv(bitflag) == False:
                adj_lpos = int(lpos)
                for tup in cigar_tups:
                    if tup[1] in ["M", "I", "N", "D"]:
                        adj_lpos += int(tup[0])

            #Now create a uniq record and if uniq store it in our set
            uniq_rec = (adj_lpos, umi, fwORrv(bitflag))
            if uniq_rec not in uniq_recs:
                #write line to file
                output.write(line)
                #add unique records to our set
                uniq_recs.add(uniq_rec)

print(uniq_recs)
'''
os.remove("temp.bam")
os.remove("sorted.sam")
'''

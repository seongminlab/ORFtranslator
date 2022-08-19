from Bio import SeqIO
from Bio.Seq import Seq
import sys,os,argparse
from itertools import product

filepath = os.path.dirname(os.path.realpath(__file__))
#print(filepath)
init_codon_list = open(filepath+"/"+ "Initiation_Codon_lists.txt","r").read().splitlines()
#print(init_codon_list)



class ORF:
    def __init__(self,sequence, start, end):
        self.seq = sequence
        self.start = start 
        self.end = end


def frameFinder(seq):

    Frame1 = [seq[i:i+3] for i in range(0, len(seq), 3)]
    Frame2 = [seq[i:i+3] for i in range(1, len(seq), 3)]
    Frame3 = [seq[i:i+3] for i in range(2, len(seq), 3)]

    Frames = [Frame1,Frame2,Frame3]

    return Frames


def findORF(Frame,init_codon,minlen,Num,init_codon_list):
    #Num: Frame+1, Frame + 2, Frame+3 : 0,1,2
    StartCodonLoci = []
    StopCodonLoci = []

    if init_codon == 0:
        for i in range(0, len(Frame)):  # Only ATG for start codon
            if Frame[i] == "ATG":
                StartCodonLoci.append(i)

    elif init_codon == 1:  # near cognate codons for start , see details: "Initiation_Codon_lists.txt" file
        for i in range(0, len(Frame)):
            if Frame[i] in init_codon_list:
                StartCodonLoci.append(i)

    elif init_codon == 2: #Any codon for start
        StartCodonLoci = [0]

    #Stop codon locus finder
    for i in range(0, len(Frame)):
        if Frame[i] == "TAA" or Frame[i] == "TAG" or Frame[i] == "TGA":
            StopCodonLoci.append(i)

#    Start Codon loci & Stop Codon loci found in frame

    combi = list(product(StartCodonLoci,StopCodonLoci))
    allow_combi = []

    for i in range(0, len(combi)):
#        print(combi)
        A = combi[i][0]
        B = combi[i][1]

        if B - A >= minlen:
            allow_combi.append(combi[i])



    not_truncated_orfs = []
    for i in range(0, len(allow_combi)):
        StartCodon = allow_combi[i][0]
        EndCodon = allow_combi[i][1]
        if not "TAA" in Frame[StartCodon:EndCodon] and not "TAG" in Frame[StartCodon:EndCodon] and not "TGA" in Frame[StartCodon:EndCodon]:
            not_truncated_orfs.append(allow_combi[i])
#    print(not_truncated_orfs)  # Remove truncated ORFs from any allow combinations


    end_codons,start_codons,uniq_end_codon_pos = [],[],[]
    if len(not_truncated_orfs) > 0:
        end_codons = list(zip(*not_truncated_orfs))[1]
        start_codons = list(zip(*not_truncated_orfs))[0]
        uniq_end_codon_pos = list(set(list(zip(*not_truncated_orfs))[1]))


    allow_orfs = []
    for i in range(0, len(uniq_end_codon_pos)):
        end_loci = uniq_end_codon_pos[i]
#        print(end_loci)
        index = list(filter(lambda x: end_codons[x] == end_loci, range(len(end_codons))))
#        print(index)
#        print(not_truncated_orfs[index[0]])
        allow_orfs.append(not_truncated_orfs[index[0]])

#    print(allow_orfs)




    orfs = []
    for i in range(0,len(allow_orfs)):
        StartCodon = allow_orfs[i][0]
        EndCodon = allow_orfs[i][1]
        if not "TAA" in Frame[StartCodon:EndCodon] and not "TAG" in Frame[StartCodon:EndCodon] and not "TGA" in Frame[StartCodon:EndCodon]:
            sequence = "".join(Frame[StartCodon:EndCodon+1])
            #print(sequence)       
            start = 3*StartCodon +1 + Num
            end = 3*(EndCodon+1) + Num
            orfs.append(ORF(sequence,start,end))
            #Start and End : 1-base rules of nucleotide positions
    return orfs




def Sequence_to_ORFs(sequence, Strand, init_codon, minlen,CodonTable,init_codon_list):
    Frames = frameFinder(sequence)
    results = []
    for i in range(0,3):
        Frame_orfs = findORF(Frames[i],init_codon,minlen,i,init_codon_list)
        Frame_number = i+1
        for orf in Frame_orfs:
            AA = str(Seq(orf.seq).translate(table=CodonTable))
            nt = str(orf.seq)
            FrameNumber = str(Frame_number)
            Strand = str(Strand)
            length = str(len(AA)-1)  # without stop * mark

            if Strand == "+":
                Start = str(orf.start)
                End = str(orf.end)

            elif Strand == "-":
                Start = str(len(sequence) - orf.start)
                End = str(len(sequence) - orf.end)

            results.append([AA,nt,Strand,FrameNumber,Start,End,length])

    return results


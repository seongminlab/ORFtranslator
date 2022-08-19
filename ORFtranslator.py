from Bio import SeqIO
from Bio.Seq import Seq
import sys,os,argparse
from itertools import combinations
from seq_utils import frameFinder,ORF,findORF,Sequence_to_ORFs

def ORFtranslator(input_fasta,minlen,init_codon,strand,CodonTable,init_codon_list):

    process_counter = 0
    combi_counter = 1
    for recod in SeqIO.parse(input_fasta,"fasta"):
        seqid = str(recod.id.split('\t')[0].split(" ")[0])
        desc = str(recod.description)
        sequence = str(recod.seq).upper().replace("U","T")
        #print(seqid,desc,sequence)

        if strand == "plus" :
            Strand = "+"

            results = Sequence_to_ORFs(sequence, Strand, init_codon, minlen,CodonTable,init_codon_list)

        elif strand == "minus":
            Strand = "-"
            sequence = str(Seq(sequence).reverse_complement())            
    
            results = Sequence_to_ORFs(sequence, Strand, init_codon, minlen,CodonTable,init_codon_list)

        elif strand == "both":
            plus_sequence = sequence
            plus_results = Sequence_to_ORFs(plus_sequence, "+", init_codon, minlen,CodonTable,init_codon_list)
            
            reverse_sequence = str(Seq(sequence).reverse_complement())
            minus_results = Sequence_to_ORFs(reverse_sequence, "-", init_codon, minlen,CodonTable,init_codon_list)

            results = plus_results + minus_results

        else:
            print("Check your parameter!! -strand  plus|minus|both")
            parser.print_help()
            sys.exit()

        
        for i in range(0,len(results)):
            Meta.write(seqid+"@orf"+str(i+1)+'\t'+'\t'.join(results[i])+'\n')
            pep.write(">"+seqid+"@orf"+str(i+1)+'\n'+str(results[i][0]).rstrip("*") +'\n')
            cds.write(">"+seqid+"@orf"+str(i+1)+'\n'+str(results[i][1]) + '\n')

            #print(seqid+"@orf"+str(i+1),'\t'.join(results[i]))

    process_counter = process_counter +1
    if process_counter == 5000:
        print(str(combi_counter*process_counter) + "sequences were processed")
        combi_counter = combi_counter+1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "ORF finder with given sequence")

    parser.add_argument("-i",metavar=":Input file" ,help="Input sequence file,FASTA format", type=str )
    parser.add_argument("-o",metavar=":Output file name",help="Output file name")
#    parser.add_argument("-type",type=str,metavar=":Input Sequence type" ,help="Sequence type, select 'DNA' or 'RNA', default is DNA",default="DNA")
    parser.add_argument("-min",metavar=":Min length",type=int,help="Minimum length of amino acid sequences, default = 10aa",default=10)
    parser.add_argument("--start",metavar=":ORF start codon to use",type=int, help="ORF start codon to use: \n 0 = ATG only, \n 1 = Alternative initiation codons (include near-cognate codons, see documents), \n 2 = any sense codon \n Default = '1'",default = 1)
    parser.add_argument("--codontable",metavar=":Genetic Code for translations",type=int, help="Select Codon Table, default = Standard Codon Table (1) \n Genetic code to use (1-31)\n see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for detail", default=1)
    parser.add_argument("--strand",metavar=":Sequence strand",help="Output ORFs on specific strand only (plus|minus|both) \n default = 'plus'",default="plus")
    
    args = parser.parse_args()

#    print(args)



    if args.i is None or args.o is None:
        parser.print_help()
        sys.exit()
    
    if args.i == args.o:
        print("\n\n*************************************************\nWarning: input file name is same to output file name!!  See usage!\n*************************************************\n\n")
        parser.print_help()
        sys.exit()

    input_fasta= args.i    
    output = args.o
    minlen = args.min
    init_codon = args.start
    CodonTable = args.codontable
    strand = args.strand

    filepath = os.path.dirname(os.path.realpath(__file__))
    init_codon_list = open(filepath+"/"+ "Initiation_Codon_lists.txt","r").read().splitlines()
    
    Meta = open(output +".meta.txt","w")
    pep = open(output+".pep.fa","w")
    cds = open(output+".cds.fa","w")


    Meta.write("ID\tPeptide\tCDS\tStrand\tFrame\tStart\tEnd\tLength_AA\n")

    ORFtranslator(input_fasta,minlen,init_codon,strand,CodonTable,init_codon_list)

    Meta.close()
    pep.close()
    cds.close()





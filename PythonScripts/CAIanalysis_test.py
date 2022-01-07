
from itertools import chain
#import Bio.Data.CodonTable as ct
from scipy.stats import gmean
from collections import Counter
import sys

inFile=open(sys.argv[1])
#test=open(sys.argv[2])
#test=test.readlines()
sequences=inFile.readlines() # list
ref={}
for i in range(0,len(sequences),2):
    ref[sequences[i].strip()[1:]]=sequences[i+1].strip()

# count the number of each codon in the sequences
print(list(ref.items())[1:5])

no_three=0
ref2={}
#print ("total transcripts are:"+str(len(ref))) #106652
for key in ref:
    seq = ref[key].upper()
    if len(seq)%3==0 and ('N' not in seq):
        codon=[]
        for i in range(0,len(seq),3):
            codon.append(seq[i:i+3])
        ref2[key]=codon
    else:
        no_three+=1

print("total transcripts that have remainder when divided by 3 or has N in sequence:"+str(no_three)) #21422

#counts have fequencey of each codon in all transcipts.
#codons = chain.from_iterable(ref2) # flat list of all codons (to be used for counting)
#counts = Counter(codons)
#print (counts) #frequency of each codon


#read a table of codon score coefficient.
csc=open(sys.argv[2])
csc_=csc.readlines()
table={}
for line in csc_[1:]:
    line=line.strip()
    table[line.split("\t")[0]]=line.split("\t")[2]

print (table)
stop_codons=['TAA', 'TAG', 'TGA']
CSC_mean={}
for key in ref2:
    cscsum=0
    for codon in ref2[key]:
        if codon not in stop_codons:
            cscsum+=float(table[codon])
    CSC_mean[key]=cscsum/len(ref2[key])
print(list(CSC_mean.items())[1:5])

outFile=open(sys.argv[3],"w")
outFile.write("ID\tMean_CSC\n")
for key in CSC_mean:
    outFile.write(key+"\t"+str(CSC_mean[key])+"\n")
outFile.close()

#print (ct.unambiguous_dna_by_id[11].forward_table)
''' register_ncbi_table(name='Bacterial, Archaeal and Plant Plastid',
 720                      alt_name=None, id=11,
 721                      table={'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
 722                             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
 723                             'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
 724                             'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
 725                             'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
 726                             'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
 727                             'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
 728                             'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
 729                             'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
 730                             'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
 731                             'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
 732                             'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
 733                             'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
 734                             'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
 735                             'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
 736                             'GGG': 'G'},
 737                      stop_codons=['TAA', 'TAG', 'TGA'],
 738                      start_codons=['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG',
 739                                    'GTG']) '''

# "if a certain codon is never used in the reference set... assign [its
# count] a value of 0.5" (page 1285)
genetic_code=1
for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table: #key
    if counts[codon] == 0:
        counts[codon] = 0.5

# determine the synonymous codons for the genetic code
def _synonymous_codons(genetic_code_dict):

    # invert the genetic code dictionary to map each amino acid to its codons
    codons_for_amino_acid = {}
    for codon, amino_acid in genetic_code_dict.items():
        codons_for_amino_acid[amino_acid] = codons_for_amino_acid.get(amino_acid, []) #Value to be returned if the key is not found
        codons_for_amino_acid[amino_acid].append(codon)

    # create dictionary of synonymous codons
    # Example: {'CTT': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'ATG': ['ATG']...}
    return {codon: codons_for_amino_acid[genetic_code_dict[codon]] for codon in genetic_code_dict.keys()}




#print(ct.unambiguous_dna_by_id.items())

_synonymous_codons = {k: _synonymous_codons(v.forward_table) for k, v in ct.unambiguous_dna_by_id.items()}
_non_synonymous_codons = {k: {codon for codon in v.keys() if len(v[codon]) == 1} for k, v in _synonymous_codons.items()}

synonymous_codons = _synonymous_codons[genetic_code]
# hold the result as it is being calulated
result = {}

 # calculate RSCU values: Relative frequency of codon to expected frequency which is mean  of all synonymous codons
for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table:
    result[codon] = counts[codon] / ((len(synonymous_codons[codon]) ** -1) * (sum((counts[_codon] for _codon in synonymous_codons[codon]))))

# calculate the weights
weights = {}
RSCUs=result #relative synonymous codon usage (RSCU) for a set of sequences
for codon in RSCUs:
        weights[codon] = RSCUs[codon] / max((RSCUs[_codon] for _codon in synonymous_codons[codon]))

print(weights) # frequency of codon RSCU/ most frequently used codon RSCU #Dictinary

# CAI calculation for each sequence #list.
# list
test_ref=[]
test_id=[]
for i in range(len(test)):
    if test[i].strip()[0]!=">":
        test_ref.append(test[i].strip())
    else:
        test_id.append(test[i].strip())




for x in range(len(test_ref)):
    sequence=test_ref[x]
    sequence=sequence.upper()
    if len(sequence)%3==0 and ('N' not in sequence):
        sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)] #list if codons in a CDS.
        sequence_weights = []
        for codon in sequence:
            if codon not in _non_synonymous_codons[genetic_code]:
                try:
                    sequence_weights.append(weights[codon])
                except KeyError:
                    # ignore stop codons
                    if codon in ct.unambiguous_dna_by_id[genetic_code].stop_codons:
                        pass
                    else:
                        raise KeyError("Bad weights dictionary passed: missing weight for codon.")

        print(str(test_id[x][1:])+'\t'+str(gmean(sequence_weights))+'\t')





'''
if sequences:
        RSCUs = RSCU(sequences, genetic_code=genetic_code)

    # determine the synonymous codons for the genetic code
    synonymous_codons = _synonymous_codons[genetic_code]

    # calculate the weights
    weights = {}
    for codon in RSCUs:
        weights[codon] = RSCUs[codon] / max((RSCUs[_codon] for _codon in synonymous_codons[codon]))

    return weights
'''

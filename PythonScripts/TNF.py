
from itertools import chain
import Bio.Data.CodonTable as ct
from scipy.stats import gmean
from collections import Counter
import numpy
import sys
inFile=open(sys.argv[1])
test=open(sys.argv[2])
outfile=open(sys.argv[3],"w")
test=test.readlines()
sequences=inFile.readlines() # list
ref=[]
id=[]
for i in range(len(sequences)):
    if sequences[i].strip()[0]!=">":
        ref.append(sequences[i].strip())
    else:
        id.append(sequences[i].strip())

# count the number of each codon in the sequences
print (len(ref))
print (len(id))
ref2=[]
no_three=0
#print ("total transcripts are:"+str(len(ref))) #106652
for seq in ref:
    seq = seq.upper()
    if len(seq)%3==0 and ('N' not in seq):
        codon=[]
        for i in range(0,len(seq),3):
            codon.append(seq[i:i+3])
        ref2.append(codon)
    else:
        no_three+=1

print("total transcripts that have remainder when divided by 3 or has N in sequence:"+str(no_three)) #21422

#counts have fequencey of each codon in all transcipts.
codons = chain.from_iterable(ref2) # flat list of all codons (to be used for counting)
counts = Counter(codons)
print (counts) #frequency of each codon

#print (ct.unambiguous_dna_by_id[11].forward_table)

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
    if RSCUs[codon]==max((RSCUs[_codon] for _codon in synonymous_codons[codon])):
        weights[codon]=1
    else:
        weights[codon]=0

print(weights) # frequency of codon RSCU/ most frequently used codon RSCU #Dictinary
'''
for codon in RSCUs:
        weights[codon] = RSCUs[codon] / max((RSCUs[_codon] for _codon in synonymous_codons[codon]))

print(weights)
'''
# CAI calculation for each sequence #list.
# list
test_ref=[]
test_id=[]
for i in range(len(test)):
    if test[i].strip()[0]!=">":
        test_ref.append(test[i].strip())
    else:
        test_id.append(test[i].strip())



'''
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

    print(str(test_id[x][1:])+'\t'+str(numpy.mean(sequence_weights))+'\t')

'''


def create_bins(lower_bound, width, quantity):
    """ create_bins returns an equal-width (distance) partitioning.
        It returns an ascending list of tuples, representing the intervals.
        A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0
        and i < quantity, satisfies the following conditions:
            (1) bins[i][0] + width == bins[i][1]
            (2) bins[i-1][0] + width == bins[i][0] and
                bins[i-1][1] + width == bins[i][1]
    """

    bins = []
    for low in range(lower_bound,
                     lower_bound + quantity * width + 1, width):
        bins.append((low, low + width))
    return bins

bins = create_bins(lower_bound=0,
                   width=1,
                   quantity=100)


def find_bin(value, bins):
    """ bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
        binning returns the smallest index i of bins so that
        bin[i][0] <= value < bin[i][1]
    """

    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1

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
            else:
                sequence_weights.append(1)

        if len(sequence_weights)>=100:

            outfile.write((str(test_id[x][1:])))
            binned_weights = {}
            for i in range(len(sequence_weights)):
                n=len(sequence_weights)
                newindex=(i*100)/float(n)
                bin_index=find_bin(newindex,bins)
                if bin_index not in binned_weights:
                    binned_weights[bin_index]=[sequence_weights[i]]
                else:
                    binned_weights[bin_index].append(sequence_weights[i])

            for key,value in binned_weights.items():
                outfile.write('\t'+str(numpy.mean(value)))

        outfile.write('\n')

outfile.close()

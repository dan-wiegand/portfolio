# Author: Daniel Wiegand
# Institution: HMS, Wyss Institute at Harvard University

from Bio.Seq import Seq
from Bio import SeqIO
import collections
import pandas
import primer3
import warnings

def GenerateMutPrimers(sequence, pos, mutation):

    position = pos - 1

    # Import Amino Acid Sequence from FASTA File
    wt_sequence = SeqIO.parse(open(sequence), 'fasta')
    sequences = collections.OrderedDict()
    for seq in wt_sequence:
        id, dna = seq.description, str(seq.seq)
        sequences[id] = dna.lower()
    active_sequence = ((list(sequences.items())[0])[1])

    # Convert Amino Acid sequence to DNA Sequence to generate primers based on E.coli codon table
    amino_acids = ["G","A","V","C","P","L","I","M","W","F","K","R","H","S","T","Y","N","Q","D","E"]
    triplets = ["GGC","GCG","GTT","TGC","CCG","CTG","ATT","ATG","TGG","TTT","AAA","CGT","CAT","AGC","ACC","TAT","AAC","CAG","GAT","GAA"]

    # Break sequence up in triplets, translate into codons, and place into ordered dictionary
    n = 3
    broken_sequence = [active_sequence[i:i+n] for i in range(0, len(active_sequence), n)]

    codons = []
    for triplet in broken_sequence:
        codon = Seq(triplet).translate(table=1,to_stop=False)
        codons.append(str(codon))

    wt_position = codons[position]

    # Generate upstream & downstream of the mutant dna primer based on requested mutagenesis position
    upstream_dna = broken_sequence[position-3] + broken_sequence[position-2] + broken_sequence[position-1]
    downstream_dna = broken_sequence[position+1] + broken_sequence[position+2] + broken_sequence[position+3]


    # Generate the mutant DNA sequence based on the requested amino acid change (triplets are hard coded from e.coli codon table

    my_mutation_dict = dict(zip(amino_acids,triplets))

    #keys,values
    # Using the codon table, generate the triplicate for the requested mutation
    mut_primer = upstream_dna + str(my_mutation_dict.get(mutation)) +  downstream_dna
    print(mut_primer)

    # Keep original reverse primer before looking into modifying it
    o_reverse_primer_list = []
    for i in range(0, 7, 1):
        o_reverse_primer_list.append(broken_sequence[position - (4 + i)])
    o_reverse_primer_list.reverse()
    reverse_original = Seq(''.join(o_reverse_primer_list)).reverse_complement()

    output = []
    x = 7
    a = 0

    while not output:
        print(a)

        reverse_primer_list = []
        for i in range(0, x, 1):
            reverse_primer_list.append(broken_sequence[position - (4 + i)])
        reverse_primer_list.reverse()
        generated_reverse_primer = Seq(''.join(reverse_primer_list)).reverse_complement()


        # perform thermodynamics analysis on forward and reverse primers
        r_tm, r_dG = primer3.calcTm(str(generated_reverse_primer).upper(), 50, 1.5, 0.5), primer3.calcHairpin(str(generated_reverse_primer).upper(), 50, 1.5, 0.5)
        f_tm, f_dG = primer3.calcTm(mut_primer.upper(), 50, 1.5, 0.5), primer3.calcHairpin(mut_primer.upper(), 50, 1.5, 0.5)

        # Calculate the absolute difference between the primer melting temperatures
        t = abs((r_tm - f_tm))
        print(t)

        # Check to see if the difference in melting temperatures if greater than a set number (n = 3)
        # If n is too large, generate a new sequence with an addition sequence to test if that is better
        # If after three attempts, this does not find something try lower the number instead
        # Else append the output as usual

        print(str(generated_reverse_primer))

        if t > 4: # Check if absolute temperature difference is greater than accepted value.
            if a == 0: # Check if we should be increasing the sequence after trying
                if x >= 10: # If we have tried 3 times to fix the difference in melting temperatures, move on
                    print("Attempted multiple times to resolve melting temperature difference by increasing sequence length.")
                    a = 1
                    x = 7
                else: # Extend the sequence by 3 to try to change melting temperature
                    x += 1
                    print("Extended by 1")
            elif a == 1: # Check if we should be decreasing the sequence
                if x <= 4: # If we have tried 3 times to fix the difference in melting temperature, move on
                    print("Attempted multiple times to resolve melting temperature difference by decreasing sequence length.")
                    x = 7
                    a = 2
                else: # Decrease the sequence by 3 to try to change melting temperature
                    x -= 1
                    print("Decreased by 1")
            elif a == 2: # If we could not fix the problem by either increasing/decreasing
                warnings.warn("Could not resolve the difference in melting temperature by changing the sequence length.")
                output.append(wt_position)
                output.append(mut_primer)
                output.append(str(reverse_original))
                output.append(t)
                output.append(f_tm)
                output.append(r_tm)
                # Need to fix this to take the lowest value
        else:
            output.append(wt_position)
            output.append(mut_primer)
            output.append(str(generated_reverse_primer))
            output.append(t)
            output.append(f_tm)
            output.append(r_tm)


    return output

sequence = "Schizo_TUTase_Codon_Optimization.fasta"
mutation_positions = [101,103,330,333,88,336,193,197,171,172,210,211,90,338,212]
mutation_positions.sort()
amino_acids = ["G","A","V","C","P","L","I","M","W","F","K","R","H","S","T","Y","N","Q","D","E"]

AA_change = []
reverse_primer = []
forward_primer = []
tm_change = []
f_temp = []
r_temp = []

for i in mutation_positions:
    for t in amino_acids:
        mutation = GenerateMutPrimers(sequence, i, t)
        AA_change.append(str("S_pombe_CID1_PuP_") + mutation[0] + str(i) + str(t))
        reverse_primer.append(mutation[1])
        forward_primer.append(mutation[2])
        tm_change.append(mutation[3])
        f_temp.append(mutation[4])
        r_temp.append(mutation[5])

mutation_data = list(zip(AA_change,reverse_primer,forward_primer,tm_change,f_temp,r_temp))

df = pandas.DataFrame(data=mutation_data, columns=['Change', 'Forward', 'Reverse','Delta_Tm', "F Temp", "R_Temp"])

print(df)

df.to_csv('PuP_Mutagenesis_Output.csv',index=False,header=True)


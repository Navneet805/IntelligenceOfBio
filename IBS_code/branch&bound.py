from Bio.Seq import Seq

text = ("GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGA"
    "TAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTC"
    "TCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTGGAAGAAATGG")

dna_seq = Seq(text)
protein_seq = str(dna_seq.translate(to_stop=False)).replace("*", "")

print("Protein Sequence:")
print(protein_seq)

AA_MW = {
'G':57, 'A':71, 'S':87, 'P':97, 'V':99,
'T':101, 'C':103, 'I':113, 'L':113, 'N':114,
'D':115, 'K':128, 'Q':128, 'E':129, 'M':131,
'H':137, 'F':147, 'R':156, 'Y':163, 'W':186
}

k = 4

def branch_and_bound(protein, k):
    max_weight = 0
    best = ""

    for i in range(len(protein)-k+1):
        weight = 0

        for j in range(k):
            weight += AA_MW[protein[i+j]]
            if weight < max_weight and j == k-1:
                break

        if weight > max_weight:
            max_weight = weight
            best = protein[i:i+k]

    return best, max_weight


print("Branch & Bound:", branch_and_bound(protein_seq, k))


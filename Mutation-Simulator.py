import random


GENETIC_CODE = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TGT': 'C', 'TGC': 'C',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'TTT': 'F', 'TTC': 'F',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAT': 'H', 'CAC': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'AAA': 'K', 'AAG': 'K',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATG': 'M',
    'AAT': 'N', 'AAC': 'N',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGA': '*'  # Stop codons
}

def create_random_gene(length=300):
    """Generates a random DNA sequence of a given length."""
    nucleotides = ['A', 'T', 'C', 'G']
    return "".join(random.choices(nucleotides, k=length))

def translate_gene(gene):
    """Translates a DNA sequence into an amino acid sequence."""
    if len(gene) % 3 != 0:
        return "Invalid gene length for translation."

    protein_sequence = []
    for i in range(0, len(gene), 3):
        codon = gene[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, 'X')  # 'X' for unknown codons
        protein_sequence.append(amino_acid)

    return "".join(protein_sequence)

def introduce_mutation(gene):
    """Introduces a single random mutation (substitution, insertion, or deletion) into the gene."""
    gene_list = list(gene)
    length = len(gene_list)
    mutation_type = random.choice(['substitution', 'insertion', 'deletion'])
    position = random.randint(0, length - 1)
    nucleotides = ['A', 'T', 'C', 'G']

    if mutation_type == 'substitution':
        original_nuc = gene_list[position]
        new_nuc = random.choice([n for n in nucleotides if n != original_nuc])
        gene_list[position] = new_nuc
        mutated_gene = "".join(gene_list)
        return mutated_gene, "substitution"

    elif mutation_type == 'insertion':
        new_nuc = random.choice(nucleotides)
        gene_list.insert(position, new_nuc)
        mutated_gene = "".join(gene_list)
        return mutated_gene, "insertion"

    elif mutation_type == 'deletion':
        del gene_list[position]
        mutated_gene = "".join(gene_list)
        return mutated_gene, "deletion"

def classify_mutation(original_gene, mutated_gene, mutation_nature):
    """Compares original and mutated proteins to classify the mutation."""
    original_protein = translate_gene(original_gene)
    mutated_protein = translate_gene(mutated_gene)

    print(f"Original DNA: {original_gene[:50]}...")
    print(f"Mutated DNA : {mutated_gene[:50]}...")
    print("-" * 50)
    print(f"Original Protein: {original_protein[:15]}...")
    print(f"Mutated Protein : {mutated_protein[:15]}...")
    print("-" * 50)
    print(f"Mutation Nature: {mutation_nature.capitalize()}")
    print("-" * 50)

    if mutation_nature in ['insertion', 'deletion'] and len(original_gene) % 3 != len(mutated_gene) % 3:
        print("Mutation Type: **Frameshift Mutation**")
        print("Explanation: An insertion or deletion shifts the reading frame, altering all downstream codons.")
    elif original_protein == mutated_protein:
        print("Mutation Type: **Silent Mutation**")
        print("Explanation: The codon change does not alter the amino acid, resulting in the same protein.")
    elif '*' in mutated_protein and mutated_protein.index('*') < len(original_protein):
        print("Mutation Type: **Nonsense Mutation**")
        print("Explanation: The codon is changed to a stop codon (*), leading to a prematurely terminated protein.")
    else:
        print("Mutation Type: **Missense Mutation**")
        print("Explanation: The codon change results in a different amino acid, altering the protein sequence.")

def main():
    """Main function to run the simulation."""
    gene_length = 300
    original_gene = create_random_gene(gene_length)
    
    
    if len(original_gene) % 3 != 0:
        original_gene = original_gene[:-(len(original_gene) % 3)]

    mutated_gene, mutation_nature = introduce_mutation(original_gene)

    if mutation_nature in ['insertion', 'deletion'] and len(mutated_gene) % 3 != 0:
        pass 

    classify_mutation(original_gene, mutated_gene, mutation_nature)

if __name__ == "__main__":
    main()


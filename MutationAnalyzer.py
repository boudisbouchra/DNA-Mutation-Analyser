class DNASequence:
    """Class representing a DNA sequence with validation and translation capabilities."""

    GENETIC_CODE = {
        "ATA": "I","ATC": "I","ATT": "I","ATG": "M","ACA": "T","ACC": "T","ACG": "T","ACT": "T","AAC": "N",
        "AAT": "N","AAA": "K","AAG": "K","AGC": "S","AGT": "S","AGA": "R","AGG": "R","CTA": "L","CTC": "L",
        "CTG": "L","CTT": "L","CCA": "P","CCC": "P","CCG": "P","CCT": "P","CAC": "H","CAT": "H","CAA": "Q",
        "CAG": "Q","CGA": "R","CGC": "R","CGG": "R","CGT": "R","GTA": "V","GTC": "V","GTG": "V","GTT": "V",
        "GCA": "A","GCC": "A","GCG": "A","GCT": "A","GAC": "D","GAT": "D","GAA": "E","GAG": "E","GGA": "G",
        "GGC": "G","GGG": "G","GGT": "G","TCA": "S","TCC": "S","TCG": "S","TCT": "S","TTC": "F","TTT": "F",
        "TTA": "L","TTG": "L","TAC": "Y","TAT": "Y","TAA": "*","TAG": "*","TGA": "*"}

    def __init__(self, sequence: str, name: str = "Sequence"):
        self.name = name
        self.sequence = sequence.upper().strip()
        self.validate_sequence()

    def validate_sequence(self):
        valid_bases = {"A", "T", "C", "G"}
        invalid_bases = set(self.sequence) - valid_bases
        if invalid_bases:
            raise ValueError(f"Invalid bases {invalid_bases} found in {self.name}.")

    def translate(self) -> str:
        protein = []
        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i : i + 3]
            amino_acid = self.GENETIC_CODE.get(codon, "?")
            if amino_acid == "*":
                protein.append("*")
                break
            protein.append(amino_acid)
        return "".join(protein)


class MutationAnalyzer:
    """Class responsible for analyzing mutations and classifying their biological impact."""

    def __init__(self, wild_type: DNASequence, mutated: DNASequence):
        self.wild_type = wild_type
        self.mutated = mutated

    def calculate_levenshtein(self) -> int:
        seq1, seq2 = self.wild_type.sequence, self.mutated.sequence
        len1, len2 = len(seq1), len(seq2)

        dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]

        for i in range(len1 + 1):
            dp[i][0] = i
        for j in range(len2 + 1):
            dp[0][j] = j

        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
                dp[i][j] = min(
                    dp[i - 1][j] + 1,  # Deletion
                    dp[i][j - 1] + 1,  # Insertion
                    dp[i - 1][j - 1] + cost,  # Substitution
                )

        return dp[len1][len2]

    def classify_impact(self) -> dict:
        wt_dna = self.wild_type.sequence
        mut_dna = self.mutated.sequence
        wt_prot = self.wild_type.translate()
        mut_prot = self.mutated.translate()
        distance = self.calculate_levenshtein()

        length_diff = abs(len(wt_dna) - len(mut_dna))

        # Check for Indels / Frameshifts
        if length_diff != 0:
            if length_diff % 3 != 0:
                classification = (
                    "FRAMESHIFT MUTATION (Indel alters reading frame)"
                )
            else:
                classification = "IN-FRAME INDEL (Insertion/Deletion of complete codon)"
        elif wt_dna == mut_dna:
            classification = "NO MUTATION (Identical sequences)"
        elif wt_prot == mut_prot:
            classification = (
                "SILENT MUTATION (DNA changed, Protein unchanged)"
            )
        elif "*" in mut_prot and mut_prot.find("*") < len(wt_prot):
            classification = (
                "NONSENSE MUTATION (Premature STOP codon introduced)"
            )
        else:
            classification = (
                "MISSENSE MUTATION (Amino acid sequence altered)"
            )

        return {
            "Levenshtein Distance": distance,
            "WT Protein": wt_prot,
            "Mutated Protein": mut_prot,
            "Classification": classification,
        }


# --- Verification Tests ---
if __name__ == "__main__":
    # Test 1: Silent Mutation (GCA -> GCC both encode Alanine)
    wt1 = DNASequence("ATGGCA", name="WT")
    mut1 = DNASequence("ATGGCC", name="Mutated")
    analyzer1 = MutationAnalyzer(wt1, mut1)

    # Test 2: Missense Mutation (GCA [Ala] -> GTA [Val])
    wt2 = DNASequence("ATGGCA", name="WT")
    mut2 = DNASequence("ATGGTA", name="Mutated")
    analyzer2 = MutationAnalyzer(wt2, mut2)

    print("Test 1 Result:", analyzer1.classify_impact()["Classification"])
    print("Test 2 Result:", analyzer2.classify_impact()["Classification"])
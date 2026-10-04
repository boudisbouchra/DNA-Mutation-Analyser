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
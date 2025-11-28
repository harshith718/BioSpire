# biospire_seq.py
# BioSpire Day-1: sequence loader, translator, basic mutation ops, and simple tests
# Requirements: python3, biopython
# Install: pip install biopython

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import json
import sys
from typing import Tuple

# ---------- Utilities ----------
def load_fasta(path: str) -> SeqRecord:
    """Load first sequence from a FASTA file."""
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError("No sequences found in FASTA: " + path)
    return records[0]

def save_fasta(seq_record: SeqRecord, path: str):
    SeqIO.write(seq_record, path, "fasta")

def clean_sequence(seq: str) -> str:
    return "".join([b for b in seq.upper() if b in "ATGC"])

# ---------- Conversions ----------
def dna_to_protein(dna_seq: str) -> str:
    return str(Seq(dna_seq).translate(to_stop=True))

def gc_content(dna_seq: str) -> float:
    dna = dna_seq.upper()
    if len(dna) == 0: return 0.0
    return 100.0 * (dna.count("G") + dna.count("C")) / len(dna)

# ---------- Mutation operators ----------
def substitute_base(seq: str, idx: int, new_base: str) -> str:
    return seq[:idx] + new_base + seq[idx+1:]

def insert_base(seq: str, idx: int, new_base: str) -> str:
    return seq[:idx] + new_base + seq[idx:]

def delete_base(seq: str, idx: int) -> str:
    return seq[:idx] + seq[idx+1:]

def random_substitution(seq: str, avoid_stop_codons: bool=True) -> Tuple[str,int,str]:
    """Randomly pick a position and substitute to a different base."""
    pos = random.randrange(len(seq))
    bases = ["A","T","G","C"]
    old = seq[pos]
    choices = [b for b in bases if b != old]
    new = random.choice(choices)
    mutated = substitute_base(seq, pos, new)
    if avoid_stop_codons:
        prot = dna_to_protein(mutated)
        if "*" in prot:
            tries = 0
            while tries < 3:
                new = random.choice(choices)
                mutated = substitute_base(seq, pos, new)
                if "*" not in dna_to_protein(mutated):
                    break
                tries += 1
    return mutated, pos, new

def random_insertion(seq: str) -> Tuple[str,int,str]:
    pos = random.randrange(len(seq)+1)
    new = random.choice(["A","T","G","C"])
    return insert_base(seq, pos, new), pos, new

def random_deletion(seq: str) -> Tuple[str,int,str]:
    if len(seq) <= 1:
        return seq, 0, ""
    pos = random.randrange(len(seq))
    return delete_base(seq, pos), pos, seq[pos]

# ---------- Small test harness ----------
def example_run():
    seq = "ATGGCTGCTGCTGAAATGGCATGCTGCTAGCTGACT"
    seq = clean_sequence(seq)
    print("Original DNA:", seq)
    print("Original protein:", dna_to_protein(seq))
    print("Original GC%:", round(gc_content(seq),2))
    s2, pos, new = random_substitution(seq)
    print(f"Substitution at {pos} -> {new}")
    print("Mutated DNA:", s2)
    print("Mutated protein:", dna_to_protein(s2))
    print("Mutated GC%:", round(gc_content(s2),2))
    s3, posi, ins = random_insertion(seq)
    print(f"Insertion at {posi} -> {ins} ; len now {len(s3)}")
    s4, posd, old = random_deletion(seq)
    print(f"Deletion at {posd} removed {old} ; len now {len(s4)}")
    rec = SeqRecord(Seq(seq), id="PETase_example", description="toy")
    save_fasta(rec, "example_petase.fasta")
    print("Saved example_petase.fasta")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        example_run()
    else:
        print("Usage: python biospire_seq.py test")

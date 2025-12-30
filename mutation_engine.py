# BioSpire: mutation utilities for evolutionary sequence simulations

import random
from typing import Tuple, Dict, List

BASES = ["A", "T", "G", "C"]

# ---------- Utilities ----------
def clean_seq(seq: str) -> str:
    """Normalize sequence by removing non-ATGC characters."""
    return "".join(c for c in seq.upper() if c in "ATGC")

def set_seed(seed: int):
    """Set random seed for reproducibility."""
    random.seed(seed)

# ---------- Mutation selection ----------
def choose_mutation_type(mutation_rates: Dict[str, float]) -> str:
    """
    Choose a mutation type based on cumulative probabilities.
    Expected keys: 'sub', 'ins', 'del'
    """
    r = random.random()
    cumulative = 0.0
    for mtype, rate in mutation_rates.items():
        cumulative += rate
        if r <= cumulative:
            return mtype
    return "sub"  # fallback

# ---------- Mutation operators ----------
def random_substitution(seq: str) -> Tuple[str, int, str]:
    if len(seq) == 0:
        return seq, -1, ""
    pos = random.randrange(len(seq))
    old = seq[pos]
    choices = [b for b in BASES if b != old]
    new = random.choice(choices)
    return seq[:pos] + new + seq[pos + 1:], pos, new

def random_insertion(seq: str) -> Tuple[str, int, str]:
    pos = random.randrange(len(seq) + 1)
    base = random.choice(BASES)
    return seq[:pos] + base + seq[pos:], pos, base

def random_deletion(seq: str) -> Tuple[str, int, str]:
    if len(seq) <= 1:
        return seq, -1, ""
    pos = random.randrange(len(seq))
    return seq[:pos] + seq[pos + 1:], pos, seq[pos]

# ---------- Mutation application ----------
def apply_mutation(seq: str, mtype: str) -> Dict:
    """Apply a single mutation of the specified type."""
    if mtype == "sub":
        new_seq, pos, base = random_substitution(seq)
        return {"seq": new_seq, "type": "sub", "pos": pos, "base": base}

    if mtype == "ins":
        new_seq, pos, base = random_insertion(seq)
        return {"seq": new_seq, "type": "ins", "pos": pos, "base": base}

    if mtype == "del":
        new_seq, pos, base = random_deletion(seq)
        return {"seq": new_seq, "type": "del", "pos": pos, "base": base}

    return {"seq": seq, "type": "none", "pos": -1, "base": ""}

def mutate_with_rate(
    seq: str,
    mutation_rates: Dict[str, float],
    mutation_prob: float
) -> Dict:
    """Apply a mutation probabilistically to a sequence."""
    if random.random() > mutation_prob:
        return {"seq": seq, "type": "none", "pos": -1, "base": ""}

    mtype = choose_mutation_type(mutation_rates)
    return apply_mutation(seq, mtype)

def batch_mutate(
    population: List[str],
    mutation_rates: Dict[str, float],
    mutation_prob: float
) -> Tuple[List[str], List[Dict]]:
    """
    Apply mutations across a population.
    Returns mutated sequences and mutation logs.
    """
    new_population = []
    logs = []

    for seq in population:
        result = mutate_with_rate(seq, mutation_rates, mutation_prob)
        new_population.append(result["seq"])
        logs.append(result)

    return new_population, logs

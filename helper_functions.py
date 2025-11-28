# mutation_engine.py
# Robust mutation utilities for BioSpire Day-2
import random
from typing import Tuple, Dict, List

BASES = ["A","T","G","C"]

def clean_seq(seq: str) -> str:
    return "".join([c for c in seq.upper() if c in "ATGC"])

def set_seed(s: int):
    random.seed(s)

def choose_mutation_type(mutation_rates: Dict[str, float]) -> str:
    r = random.random()
    cum = 0.0
    for k,v in mutation_rates.items():
        cum += v
        if r <= cum:
            return k
    return "sub"

def random_substitution(seq: str) -> Tuple[str,int,str]:
    if len(seq) == 0:
        return seq, -1, ""
    pos = random.randrange(len(seq))
    old = seq[pos]
    choices = [b for b in BASES if b != old]
    new = random.choice(choices)
    new_seq = seq[:pos] + new + seq[pos+1:]
    return new_seq, pos, new

def random_insertion(seq: str) -> Tuple[str,int,str]:
    pos = random.randrange(len(seq)+1)
    base = random.choice(BASES)
    return seq[:pos] + base + seq[pos:], pos, base

def random_deletion(seq: str) -> Tuple[str,int,str]:
    if len(seq) <= 1:
        return seq, -1, ""
    pos = random.randrange(len(seq))
    return seq[:pos] + seq[pos+1:], pos, seq[pos]

def apply_mutation(seq: str, mtype: str) -> Dict:
    if mtype == "sub":
        new_seq, pos, base = random_substitution(seq)
        return {"seq": new_seq, "type":"sub", "pos":pos, "base":base}
    if mtype == "ins":
        new_seq, pos, base = random_insertion(seq)
        return {"seq": new_seq, "type":"ins", "pos":pos, "base":base}
    if mtype == "del":
        new_seq, pos, base = random_deletion(seq)
        return {"seq": new_seq, "type":"del", "pos":pos, "base":base}
    return {"seq": seq, "type":"none", "pos":-1, "base":""}

def mutate_with_rate(seq: str, mutation_rates: Dict[str,float], mutation_prob: float) -> Dict:
    if random.random() > mutation_prob:
        return {"seq": seq, "type":"none", "pos":-1, "base":""}
    mtype = choose_mutation_type(mutation_rates)
    return apply_mutation(seq, mtype)

def batch_mutate(pop: List[str], mutation_rates: Dict[str,float], mutation_prob: float) -> Tuple[List[str], List[Dict]]:
    new_pop = []
    logs = []
    for s in pop:
        res = mutate_with_rate(s, mutation_rates, mutation_prob)
        new_pop.append(res["seq"])
        logs.append(res)
    return new_pop, logs

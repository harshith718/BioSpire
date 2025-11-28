# evolution_engine.py
import json, statistics, os
from typing import List, Dict, Any
from mutation_engine import batch_mutate, set_seed, clean_seq
from biospire_seq import dna_to_protein

# Placeholder scoring (we'll replace with scoring module later)
def sequence_score(seq: str) -> float:
    if len(seq) == 0:
        return 0.0
    gc = (seq.count("G") + seq.count("C")) / max(1, len(seq))
    score = 1.0 - abs(gc - 0.5)  # closer to 0.5 better
    length_penalty = max(0.0, 1.0 - abs(len(seq) - 300) / 1000.0)
    return max(0.0, score * length_penalty)

class Population:
    def __init__(self, sequences: List[str]):
        self.sequences = [clean_seq(s) for s in sequences]

    def evaluate(self) -> List[float]:
        return [sequence_score(s) for s in self.sequences]

    def select_top_k(self, k:int) -> List[str]:
        scored = list(zip(self.sequences, self.evaluate()))
        scored.sort(key=lambda x: x[1], reverse=True)
        return [s for s,_ in scored[:k]]

def run_evolution(initial_seqs: List[str],
                  generations: int = 50,
                  pop_size: int = 50,
                  mutation_rates: Dict[str,float] = None,
                  mutation_prob: float = 0.2,
                  selection_fraction: float = 0.5,
                  seed: int = 42,
                  out_dir: str = "../simulations",
                  log_path: str = "../logs/evolution_log.json") -> Dict[str,Any]:
    if mutation_rates is None:
        mutation_rates = {"sub":0.7, "ins":0.15, "del":0.15}
    set_seed(seed)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    pop = (initial_seqs * ((pop_size // len(initial_seqs))+1))[:pop_size]
    population = Population(pop)
    history = []
    for gen in range(generations):
        new_seqs, mut_logs = batch_mutate(population.sequences, mutation_rates, mutation_prob)
        population.sequences = new_seqs
        scores = population.evaluate()
        avg_score = statistics.mean(scores)
        max_score = max(scores)
        best_idx = scores.index(max_score)
        best_seq = population.sequences[best_idx]
        # selection
        survivors = population.select_top_k(max(1, int(selection_fraction * pop_size)))
        # refill
        new_pop = survivors.copy()
        while len(new_pop) < pop_size:
            new_pop.append(survivors[len(new_pop) % len(survivors)])
        population.sequences = new_pop
        history.append({
            "generation": gen,
            "avg_score": avg_score,
            "max_score": max_score,
            "best_seq": best_seq,
            "mutations_sample": mut_logs[:5]
        })
    with open(log_path, "w") as f:
        json.dump({"history": history}, f, indent=2)
    return {"generations": generations, "final_best": history[-1]["best_seq"], "history": history}

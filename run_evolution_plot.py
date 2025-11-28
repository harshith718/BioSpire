# run_evolution_plot.py
import os, argparse, json
from evolution_engine import run_evolution
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--generations", type=int, default=40)
parser.add_argument("--pop", type=int, default=40)
parser.add_argument("--mutprob", type=float, default=0.25)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

OUT_GRAPH_DIR = "../graphs"
os.makedirs(OUT_GRAPH_DIR, exist_ok=True)

initial = ["ATGGCTGCTGCTGAAATGGCATGCTGCTAGCTGACT"]  # placeholder; replace with real FASTA later
res = run_evolution(initial_seqs=initial,
                    generations=args.generations,
                    pop_size=args.pop,
                    mutation_prob=args.mutprob,
                    seed=args.seed,
                    out_dir="../simulations",
                    log_path="../logs/evolution_log.json")

history = res["history"]
gens = [h["generation"] for h in history]
avg = [h["avg_score"] for h in history]
mx = [h["max_score"] for h in history]

plt.figure(figsize=(7,4))
plt.plot(gens, avg, label="avg_score")
plt.plot(gens, mx, label="max_score")
plt.xlabel("Generation")
plt.ylabel("Score")
plt.title("BioSpire: Fitness over Generations")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_GRAPH_DIR, "fitness_curve.png"))
plt.close()

# quick mutation frequency heatmap (positions vs generations)
max_len = max(len(s) for s in history[-1]["best_seq"])
pos_counts = np.zeros((len(history), max_len))
for i,h in enumerate(history):
    muts = h.get("mutations_sample", [])
    for m in muts:
        p = m.get("pos", -1)
        if p >= 0 and p < max_len:
            pos_counts[i,p] += 1
plt.figure(figsize=(8,5))
plt.imshow(pos_counts.T, aspect='auto', origin='lower')
plt.xlabel("Generation")
plt.ylabel("Position")
plt.title("Mutation heatmap (sampled mutations)")
plt.colorbar(label="mutation count")
plt.tight_layout()
plt.savefig(os.path.join(OUT_GRAPH_DIR, "mutation_heatmap.png"))
plt.close()

# Save final best sequence to file for portfolio
with open(os.path.join(OUT_GRAPH_DIR, "final_best_sequence.txt"), "w") as f:
    f.write(res["final_best"])

print("Done. Graphs saved to:", os.path.abspath(OUT_GRAPH_DIR))
print("Log saved to ../logs/evolution_log.json")

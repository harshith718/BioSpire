# BioSpire — Evolution Simulation Engine  

BioSpire is a compact evolutionary simulation engine designed to explore how sequences mutate, adapt, and improve over generations.  
It models mutation rates, selection pressure, evolutionary fitness scoring, and generates useful visualizations such as heatmaps and fitness curves.

This project is part of a 6-project evolutionary research portfolio (BioSpire, EON, GeneFlux, EcoLens, BioGraph, Micro-ECs).

---

## 🔥 Features

- Simulates sequence evolution under mutation and selection  
- Tracks best-performing sequences across generations  
- Generates fitness curves and mutation heatmaps  
- Visualizes evolutionary improvement over time  
- Modular Python code structure (mutation engine, evolution engine, plotting)

---

## 📁 Project Structure

Your current GitHub project places all files in the `code/` folder:

code/
biospire_seq.py
evolution_engine.py
mutation_engine.py
helper_functions.py
run_evolution_plot.py
example_petase.fasta
evolution_log.json

graphs/
fitness_curve.png
mutation_heatmap.png
final_best_sequence.txt

logs/
evolution_log.json


---

## ▶️ How to Run

### **1. Install Python (3.9 or above recommended)**  
Check if installed:
```bash
python --version
2. Install Dependencies

Inside the project folder run:

pip install numpy matplotlib

3. Run the full simulation
python code/run_evolution_plot.py


This will:

Run evolution

Generate graphs

Save the best sequence

Save logs

Outputs appear in:

graphs/

logs/
Example Output

fitness_curve.png → Fitness score progression

mutation_heatmap.png → Mutation distribution heatmap

final_best_sequence.txt → Best evolved sequence

evolution_log.json → Detailed log

# Cuckoo Catfish Optimizer (CCO) – Algorithm Architecture

This document describes the reference implementation architecture of the Cuckoo Catfish Optimizer.

---

## 1. High-Level Execution Flow

1. Initialize population
2. Evaluate fitness
3. Repeat until termination:
   - Generate parasitic candidates
   - Perform host replacement mechanism
   - Apply adaptive exploration/exploitation step
   - Evaluate updated solutions
   - Update global best
4. Return best solution

---

## 2. Main Phases

### Phase 1 – Initialization

- Randomly generate N candidate solutions within bounds.
- Evaluate objective function.
- Sort population based on fitness.
- Store global best.

**Data structures:**
- `Population (N × D)`
- `Fitness (N × 1)`
- `BestPosition (1 × D)`
- `BestFitness (scalar)`

---

### Phase 2 – Parasitic Candidate Generation

For each individual:

- Generate a new candidate using:
  - Random perturbation
  - Directed movement toward better individuals
- Control step size via adaptive parameter.

Mathematical form (conceptual):

X_new = X_i + α * rand * (X_best − X_i) + β * random_walk

Where:
- α controls exploitation
- β controls exploration

---

### Phase 3 – Host Replacement Mechanism

- Compare newly generated candidate with a randomly selected host.
- If the candidate has better fitness:
  - Replace host.
- Otherwise:
  - Host survives.

This ensures competitive survival similar to brood parasitism.

---

### Phase 4 – Adaptive Parameter Control

Exploration parameter decreases over iterations:

α(t) = α0 × (1 − t / T)

This allows:
- Early iterations → exploration
- Later iterations → exploitation

---

### Phase 5 – Boundary Control

After each update:

- Clamp solutions within lower and upper bounds.
- Ensure feasibility.

---

## 3. Core Mechanisms

- Population-based search
- Parasitic replacement strategy
- Fitness-driven selection
- Adaptive coefficient decay
- Global best memory

---

## 4. Data Structures Used

| Structure | Shape | Purpose |
|------------|--------|----------|
| Population | N × D | Candidate solutions |
| Fitness | N × 1 | Objective values |
| BestPosition | 1 × D | Global best |
| BestFitness | Scalar | Best objective value |

---

## 5. Key Functions (Reference Implementation)

- `InitializePopulation(N, D, LB, UB)`
- `EvaluatePopulation(pop, objFun)`
- `GenerateCandidate(Xi, Xbest, parameters)`
- `ReplacementStrategy(pop, fitness)`
- `UpdateBest(pop, fitness)`
- `BoundaryControl(X, LB, UB)`

---

## 6. Logical Modules

1. Initialization module
2. Candidate generation module
3. Competitive replacement module
4. Adaptive control module
5. Termination and reporting module

---

## 7. Termination Criteria

- Maximum number of iterations
- Maximum function evaluations
- Convergence tolerance (optional)

---

This architecture ensures a modular, reproducible, and extensible implementation of the Cuckoo Catfish Optimizer suitable for integration into benchmarking frameworks.
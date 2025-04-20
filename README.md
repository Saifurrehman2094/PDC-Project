# PDC-Project


# Parallel SSSP Update Algorithm – Phase 1

## Group Members:
- Muhammad Owais -22i-0945
- Ali Haider -22i-0936
- Saif ur Rehman -22i-0933

## Selected Paper:
**Title:** A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks  
**Authors:** Arindam Khanda, Sriram Srinivasan, Sanjukta Bhowmick, Boyana Norris, Sajal K. Das  
**Link:** [https://drive.google.com/file/d/1Cj7u6bLfbwfjSwtZhfDwv4V5wIuiGB9y/view?usp=sharing](https://drive.google.com/file/d/1Cj7u6bLfbwfjSwtZhfDwv4V5wIuiGB9y/view?usp=sharing)

## Project Summary:
This project explores a parallel framework for efficiently updating Single-Source Shortest Paths (SSSP) in large-scale dynamic graphs. Rather than recomputing paths from scratch, the algorithm identifies affected subgraphs and selectively updates them. The framework is scalable and deployable across both CPU and GPU architectures.

## Phase 1 Deliverables:
- 📄 Literature Review and Paper Presentation (included in `Phase1_Presentation/`)
- 📑 Summary of the selected paper's contributions
- ⚙️ Proposed parallelization strategy

## Parallelization Strategy:
We plan to implement the algorithm using the following parallel computing tools:

### 🧩 MPI – Inter-node Communication
MPI will be used to distribute graph partitions across different nodes in a cluster, enabling concurrent processing and efficient coordination of subgraph updates.

### 🔄 OpenMP / OpenCL – Intra-node Parallelism
- **OpenMP** will handle shared memory parallelism within each node for CPU-based execution.
- **OpenCL** may be explored for GPU acceleration depending on hardware availability and scalability.

### 🧠 METIS – Graph Partitioning
METIS will be used to partition large graphs into smaller, balanced chunks for distribution via MPI. It minimizes edge cuts and ensures efficient communication between partitions.

## Repository Structure:
```
root/
│
├── Phase1_Presentation/
│   ├── SSSP_Parallel_Algorithm_Presentation.pptx
│   └── selected_paper.pdf
│
├── docs/
│   └── implementation_plan.md  # Optional for planning Phase 2
│
├── datasets/
│   └── [To be added in Phase 2]
│
├── Phase2_Code/
│   └── [Implementation will go here]
│
├── results/
│   └── [Performance graphs and metrics]
│
└── README.md
```

## Project Timeline:
- 📌 **Phase 1 – Literature Review & Presentation:** April 20
- 🛠️ **Phase 2 – Implementation & Demo:** May 4
- 🔄 Weekly GitHub commits and progress tracking required

## License:
This project is for academic purposes under FAST NUCES. No license required.

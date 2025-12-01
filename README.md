dan-wiegand — Computational Biology & Bioinformatics Portfolio

This repository contains selected Python and R tools developed during my work at the Wyss Institute at Harvard University, Harvard Medical School, and during advanced academic training in chemical and biological engineering.

These projects demonstrate end-to-end computational workflows for oligonucleotide library design, NGS error profiling, primer design, thermodynamic optimization, kinetic modeling, and algorithmic sequence engineering. All scripts are hand-written and reflect real research applications.

Project Index

001_OligoForge
File: 001_OligoForge
Automated Multiplex Oligonucleotide Library Constructor
Description:
A high-complexity DNA library builder used for multiplex assembly workflows.
Includes barcode-compatibility screening (Hamming distance and melting-temperature checks), homology-arm insertion, restriction-site placement, sequence fragmentation using overlapping blocks, and FASTA export for OLS array synthesis.
Applications: combinatorial DNA design, synthetic biology, high-plex library construction.

002_KineTraceR
File: 002_KineTraceR
Fluorescence Kinetics and Reaction-Rate Analyzer
Description:
R script for extracting reaction rates from fluorescence time-course data generated during in vitro protein synthesis assays.
Includes spline-based numerical differentiation, replicate averaging (mean and SEM), and structured data output.
Applications: enzymatic kinetics and assay development.

003_ReactSim
File: 003_ReactSim
Chemical Reactor ODE Simulation Script
Description:
Python model simulating cellulose conversion in a chemical reactor.
Implements multi-species reaction kinetics, SciPy ODE integration, concentration calculations, and reaction-rate evaluation.
Applications: chemical reaction engineering and numerical modeling.

004_OligoQC
File: 004_OligoQC
NGS-Based Oligonucleotide Library Error Profiling Pipeline
Description:
Comprehensive R pipeline for evaluating sequence fidelity in a 12K CustomArray oligonucleotide library sequenced on an Illumina platform.
Includes BAM parsing, mutation and indel quantification, GC-content analysis, spatial error mapping across array quadrants, permutation testing, chi-square modeling, and QC visualization (heatmaps, boxplots, distributions).
Applications: sequencing QC, array validation, high-throughput DNA synthesis evaluation.

005_MutPrimerDesigner
File: 005_MutPrimerDesigner
Automated Thermodynamics-Guided Mutagenesis Primer Designer
Description:
Python tool that generates forward and reverse primers for saturation mutagenesis across multiple amino-acid positions.
Includes codon substitution logic, primer Tm matching, dynamic primer-length adjustment, thermodynamic evaluation using primer3, and export of final primer sets to CSV.
Applications: protein engineering and enzyme mutagenesis.

Skills Demonstrated

Programming:
Python (Biopython, primer3, SciPy, numpy, collections)
R (Rsamtools, Biostrings, statistical modeling, visualization)
GitHub version control

Bioinformatics:
FASTA and BAM parsing
NGS quality control and error modeling
Barcode and primer design
Oligonucleotide library engineering

Computational and Analytical Skills:
Thermodynamic modeling
ODE solving and numerical simulation
Statistical inference and permutation testing
High-throughput data visualization

About This Repository

These tools were developed to support research involving high-plex fluorescence sequencing (FISSEQ), synthetic DNA library construction, NGS-based QC analytics, mutagenesis design, and biochemical reaction modeling.
They reflect extensive experience with scientific programming, custom algorithm development, and applied computational biology.


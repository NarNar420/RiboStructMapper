# RiboPrint Overview

## What is RiboPrint?

RiboPrint is a bioinformatics tool that bridges the gap between **1D functional genomic data** and **3D structural biology**. It takes linear ribosome density data (from Ribosome Profiling/Ribo-Seq experiments) and maps those density scores directly onto the 3D atomic coordinates of a protein structure.

By doing this, researchers can visually inspect exactly *where* on the folded protein structure ribosomes were actively translating, stalling, or slowing down during the experiment. 

## How Does it Work?

The pipeline operates in five core steps:

1. **Genomic Parsing**: Extracts the exact Coding Sequence (CDS) for a gene using a FASTA sequence and GTF annotation file.
2. **Density Extraction**: Pulls the linear ribosome density scores from a `bedGraph` file that match the extracted CDS coordinates.
3. **Sequence Alignment**: Translates the genomic DNA/RNA into a theoretical amino acid sequence and performs a global alignment against the actual sequences extracted from the physical `ATOM` records of the provided PDB file. This handles missing residues or gaps in the 3D structure.
4. **Offset & Aggregation**: 
   - Applies nucleotide-level "offsets" (e.g., shifting -12 or -15 nucleotides) to account for where the ribosome active site is relative to the protected mRNA fragment.
   - Aggregates the density of nucleotide triplets (codons) into a single score per amino acid (using mean, max, sum, or median).
5. **B-Factor Injection**: Writes the final amino acid scores into the **B-factor column** of the PDB file. 

## Why B-Factors?

In standard PDB files, the B-factor (temperature factor) column usually indicates structural uncertainty or atomic mobility. By intentionally overwriting this column with our ribosome density scores, the resulting PDB file can be opened in standard molecular visualization software (like PyMOL or ChimeraX). Users can then easily color the protein surface or backbone by B-factor to instantly see a 3D heatmap of translation density.

## Usage Environments

- **Web UI**: A user-friendly dark-mode web application for drag-and-drop processing.
- **REST API**: A programmatic interface (`/submit_job`, `/status`, `/download`) for integrating the tool into larger automated bioinformatic workflows.
- **CLI / Python Module**: Core functions can be run locally or integrated into custom Python scripts.

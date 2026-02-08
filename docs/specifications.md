This is a comprehensive architectural blueprint for your web application. This plan is designed to serve as the master specification document for development.

### **Project Title:** Ribo-Struct Mapper

**Objective:** A web-based pipeline to visualize ribosome density on 3D protein structures by translating genomic data and injecting density scores into the B-factor field of PDB files.

---

### **1. High-Level System Architecture**

The system will follow a classic **Client-Server** model with a stateless backend to handle the bioinformatics processing.

* **Frontend (Client):** A user-friendly interface for file uploads, parameter configuration, and job status visualization.
* **Backend (Server):** An API-driven server that accepts files, performs the heavy logical processing (alignment, mapping, injection), and serves the processed files back.
* **Processing Engine:** A modular logic core (likely Python-based due to robust PDB parsing libraries like Biopython) that handles the specific scientific calculations.

---

### **2. Logic Pipeline & Module Breakdown**

This is the core "Brain" of the app. The backend should be divided into these specific functional modules:

#### **Module A: The Input Parser**

* **Function 1: `parse_genomic_data(fasta, gtf, gene_id)**`
* **Logic:** parses the GTF/GFF to find the Coding Sequence (CDS) coordinates for the specific gene. It extracts the corresponding nucleotide sequence from the FASTA file.
* **Output:** A raw nucleotide string (DNA/RNA) and a coordinate map (Genomic Coordinate  Nucleotide Index).


* **Function 2: `parse_ribo_density(density_file, coordinates)**`
* **Logic:** Reads the bedGraph (or bigWig) file. It filters for scores only falling within the CDS coordinates found in Function 1.
* **Output:** A vector/array of raw density scores matching the nucleotide sequence length.


* **Function 3: `parse_pdb_structure(pdb_file)**`
* **Logic:** Parses the PDB to extract the amino acid (AA) sequence specifically from the `ATOM` records (the physical coordinates), *not* just the header. It must account for Chain IDs.
* **Output:** An object containing the PDB AA sequence, residue numbers, and chain identifiers.



#### **Module B: The Translator & Aligner (The Critical Step)**

* **Function 4: `translate_sequence(nucleotide_seq)**`
* **Logic:** Standard translation (genetic code table).
* **Output:** A "Theoretical" Amino Acid sequence derived from the genome.


* **Function 5: `generate_alignment_map(theoretical_aa, pdb_aa)**`
* **Logic:** Performs a pairwise sequence alignment (e.g., Needleman-Wunsch or Smith-Waterman).
* **Crucial Handling:**
* **Gaps:** If the PDB has missing residues (disordered regions), the alignment must introduce gaps in the PDB sequence so the genomic scores don't "slide" onto the wrong atoms.
* **Mismatches:** Allow for slight variations (e.g., mutations/SNPs) but flag high discordance.


* **Output:** A **Residue Mapping Index**: A dictionary or array linking `Theoretical_AA_Index`  `PDB_Residue_ID`.



#### **Module C: The Data Processor**

* **Function 6: `aggregate_codons(density_vector, method)**`
* **Logic:** Iterates through the raw density vector in steps of 3 (codons).
* **User Choice:** Applies the aggregation method (Mean, Max, Sum, or Median) to generate a single score per amino acid.
* **Output:** A `Base_Score_Vector` (length = length of protein).


* **Function 7: `apply_offset(base_score_vector, offset_value)**`
* **Logic:** Shifts the vector indices by the `offset_value`.
* *Example:* If offset is -15 (5 residues), the score typically at residue 50 is moved to residue 45.
* *Boundary Handling:* Discards scores that shift "off" the ends of the protein.


* **Output:** An `Offset_Score_Vector`.



#### **Module D: The PDB Injector**

* **Function 8: `inject_bfactors(pdb_structure, mapping_index, offset_score_vector)**`
* **Logic:**
1. Iterates through the `ATOM` records of the PDB.
2. For each atom, checks its Residue ID.
3. Look up the Residue ID in the `Mapping Index` to find the corresponding `Theoretical_AA_Index`.
4. Retrieves the score from the `Offset_Score_Vector`.
5. Overwrites the value in the **B-factor** column (columns 61-66 in standard PDB format) with the retrieved score.


* **Output:** A new `.pdb` file instance.



---

### **3. Infrastructure & API Specifications**

#### **API Endpoints**

1. **`POST /upload_job`**
* **Body:** Multipart form data containing:
* `.pdb` file
* Sequence files (FASTA + GTF)
* Density file (bedGraph)
* JSON Configuration: `{"offsets": [0, -40, -90], "aggregation": "mean", "chain": "A"}`


* **Response:** `job_id` (UUID).


2. **`GET /status/{job_id}`**
* **Response:** `{"status": "processing" | "completed" | "failed", "progress": "60%"}`


3. **`GET /download/{job_id}`**
* **Response:** A ZIP file containing all generated PDBs (e.g., `protein_offset_0.pdb`, `protein_offset_-40.pdb`).



#### **Data Storage Strategy**

* **Temporary Storage:** Files should be stored in a temporary directory hashed by `job_id`.
* **Cleanup:** Automated cron job to delete files older than 24 hours to save server space (as density files can be large).

---

### **4. User Interface (Frontend) Layout**

**Section 1: The Input Deck**

* **Structure Upload:** Drag & drop zone for `.pdb`.
* **Genomic Context:** Drag & drop for `.fasta` and `.gtf`/`.gff`.
* **Experimental Data:** Drag & drop for ribosome density file.

**Section 2: Configuration Panel**

* **Aggregation Method:** Dropdown menu (Average, Sum, Max).
* **Offsets:** An input field to accept comma-separated integers (e.g., `0, -30, -60`).
* **Chain Selection:** (Optional) Input for specific PDB chain ID (default: "A" or "All").

**Section 3: Action & Output**

* **"Run Alignment" Button:** First, runs the alignment and shows the user a "QC Score" (e.g., "98% Sequence Identity detected between Genomic and PDB").
* **"Generate Models" Button:** Executes the B-factor injection.
* **Download:** A "Download All" button appearing upon completion.

---

### **5. Edge Cases & Error Handling to Plan For**

1. **Coordinate Mismatch:** The GTF coordinates don't match the FASTA headers.
* *Solution:* Validate headers before processing.


2. **The "Missing Density" Problem:** The PDB file is missing a loop (residues 50-60 are absent in the structure but present in the genome).
* *Solution:* The alignment step must detect this. The density scores for genomic residues 50-60 are simply ignored/discarded because there are no atoms to write them to.


3. **Multiple Transcripts:** The GTF contains multiple isoforms for the gene.
* *Solution:* The app should either ask the user to specify a `transcript_id` or default to the longest CDS (Canonical).



### **6. Recommended Tech Stack**

* **Backend:** Python (FastAPI or Flask).
* *Why:* Unbeatable libraries for this specific task (`Biopython` for PDB/SeqIO, `Pandas` for BedGraph processing, `Bio.Align` for the sequence alignment).


* **Frontend:** React or simple HTML/JS (Bootstrap).
* **Deployment:** Docker container (for easy portability between your dev environment and any server).

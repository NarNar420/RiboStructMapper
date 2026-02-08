# 🚀 Quick Start Guide - Web Interface

## How to Use RiboStructMapper (5 Easy Steps)

### 1. Start the Server
```bash
# If using Docker (recommended)
docker-compose up -d

# Or manually
uvicorn ribostruct.web.server:app --reload
```

### 2. Open Your Browser
Navigate to: **http://localhost:8000**

### 3. Prepare Your Files
You need **4 files**:
- **PDB file** (`.pdb`) - Your protein structure
- **FASTA file** (`.fasta`) - Genomic sequence
- **GTF file** (`.gtf`) - Gene annotation
- **bedGraph file** (`.bedgraph`) - Ribosome density data

### 4. Upload and Process
1. **Upload files** - Drag & drop or click to select each file
2. **Enter offsets** - Type comma-separated values (e.g., `0,-12,-15`)
3. **Click "Process Files"** - Wait for processing to complete
4. **Watch status** - Progress updates every 2 seconds

### 5. Download Results
When status shows **"Completed"**:
- Click **"Download Results"**
- Get a ZIP file with all processed PDB files
- Open in PyMOL, ChimeraX, or any PDB viewer

## 📊 What the Output Contains

The ZIP file includes one PDB file per offset:
- `output_offset_0.pdb` - No offset (raw density)
- `output_offset_-12.pdb` - 12 nucleotide shift
- `output_offset_-15.pdb` - 15 nucleotide shift (if specified)

Each PDB file has **ribosome density scores in the B-factor column** - ready for visualization!

## 💡 Tips

- **Offsets**: Common values are `0, -12, -15` for P-site mapping
- **Processing time**: Usually takes 5-30 seconds depending on file size
- **Status**: "completed" = ready to download, "failed" = check your files
- **Results**: Automatically deleted after 24 hours to save space

## 🆘 Troubleshooting

**Upload fails?**  
→ Check file formats match requirements (`.pdb`, `.fasta`, `.gtf`, `.bedgraph`)

**Processing fails?**  
→ Ensure gene ID exists in GTF and chromosome names match between files

**Download not working?**  
→ Wait for status to show "completed" before downloading

---

That's it! Your density-mapped PDB files are ready for visualization. 🎉

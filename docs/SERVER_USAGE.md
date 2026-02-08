# Server Usage Instructions

## Starting the Server

To start the FastAPI server with auto-reload:

```bash
source .venv/bin/activate
uvicorn server:app --reload
```

The server will start at `http://127.0.0.1:8000`

## Using the Swagger UI

1. Open your browser and navigate to: `http://127.0.0.1:8000/docs`

2. You'll see the interactive API documentation (Swagger UI)

3. Click on `POST /submit_job` to expand it

4. Click the "Try it out" button

5. Upload your files:
   - **pdb_file**: Click "Choose File" and select your PDB file
   - **fasta_file**: Select your genomic FASTA file  
   - **gtf_file**: Select your GTF/GFF annotation file
   - **density_file**: Select your bedGraph density file
   - **offsets**: Enter comma-separated values (e.g., `0,-10,-15`)

6. Click "Execute"

7. The response will contain:
   - `job_id`: Unique identifier for your job
   - `status`: "queued" - job is queued for processing
   - `message`: Confirmation message

8. Check job progress:
   - Navigate to `jobs/{job_id}/` folder
   - Check `status.txt` file:
     - `uploaded` → `processing` → `completed` (or `failed`)
   - When completed, find output PDB files: `output_offset_0.pdb`, etc.

## Testing with Mock Data

To test with the provided mock data:

```bash
# Make sure server is running in one terminal
uvicorn server:app --reload

# In another terminal, run the test script
python test_server_manual.py
```

This will submit a job and monitor its progress automatically.

## API Endpoints

- `GET /` - API information
- `GET /health` - Health check
- `POST /submit_job` - Submit a new processing job
- `GET /docs` - Interactive API documentation (Swagger UI)
- `GET /redoc` - Alternative documentation

## Job Directory Structure

```
jobs/
└── {job_id}/
    ├── input.pdb          # Uploaded PDB file
    ├── input.fasta        # Uploaded FASTA file
    ├── input.gtf          # Uploaded GTF file
    ├── input.bedgraph     # Uploaded bedGraph file
    ├── params.txt         # Job parameters
    ├── status.txt         # Current status
    ├── output_offset_0.pdb      # Generated output (when complete)
    ├── output_offset_-10.pdb    # Generated output (when complete)
    ├── output_files.txt   # Log of generated files
    └── error.txt          # Error details (if failed)
```

## Status Values

- `uploaded`: Files uploaded, awaiting processing
- `queued`: Job queued for background processing
- `processing`: Currently being processed
- `completed`: Processing finished successfully
- `failed`: Processing encountered an error

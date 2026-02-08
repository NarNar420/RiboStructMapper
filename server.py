"""
RiboStructMapper FastAPI Server

This module provides the web API for the RiboStructMapper pipeline.
It handles file uploads, job management, and serves processed results.
"""

import os
import uuid
import shutil
import traceback
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, UploadFile, File, Form, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse
from pydantic import BaseModel

# Import the pipeline function
from main_cli import run_pipeline


# Initialize FastAPI application
app = FastAPI(
    title="RiboStructMapper API",
    description="Web API for mapping ribosome density onto PDB structures",
    version="1.0.0"
)

# Configuration
JOBS_DIR = Path("jobs")
JOBS_DIR.mkdir(exist_ok=True)


# ============================================================================
# Models
# ============================================================================

class JobResponse(BaseModel):
    """Response model for job submission."""
    job_id: str
    status: str
    message: Optional[str] = None


# ============================================================================
# Utility Functions
# ============================================================================

def create_job_id() -> str:
    """
    Generate a unique job ID using UUID4.
    
    Returns:
        A unique UUID string (e.g., "a3b8d1b6-0b3b-4b1a-9c1a-1a2b3c4d5e6f")
    """
    return str(uuid.uuid4())


def save_upload_file(upload_file: UploadFile, destination: Path) -> None:
    """
    Save an uploaded file to the specified destination.
    
    Args:
        upload_file: The uploaded file from FastAPI
        destination: Path where the file should be saved
    """
    with open(destination, "wb") as f:
        shutil.copyfileobj(upload_file.file, f)


def process_job(job_id: str) -> None:
    """
    Process a job in the background.
    
    This function loads the uploaded files, parses parameters, runs the
    pipeline, and updates the job status.
    
    Args:
        job_id: The unique job identifier
    """
    job_dir = JOBS_DIR / job_id
    status_file = job_dir / "status.txt"
    
    try:
        # Update status to processing
        with open(status_file, "w") as f:
            f.write("processing\n")
        
        # Define file paths
        pdb_path = str(job_dir / "input.pdb")
        fasta_path = str(job_dir / "input.fasta")
        gtf_path = str(job_dir / "input.gtf")
        bedgraph_path = str(job_dir / "input.bedgraph")
        
        # Parse parameters
        params_file = job_dir / "params.txt"
        with open(params_file, "r") as f:
            params_content = f.read()
        
        # Extract offsets from params.txt
        # Format: "offsets=0,-10,-15"
        offsets_str = ""
        for line in params_content.split("\n"):
            if line.startswith("offsets="):
                offsets_str = line.split("=", 1)[1].strip()
                break
        
        # Parse offsets into a list of integers
        offsets = [int(x.strip()) for x in offsets_str.split(",") if x.strip()]
        
        # Run the pipeline
        output_files = run_pipeline(
            pdb_path=pdb_path,
            fasta_path=fasta_path,
            gtf_path=gtf_path,
            bedgraph_path=bedgraph_path,
            offsets=offsets,
            aggregation_method='mean',
            output_dir=str(job_dir)
        )
        
        # Update status to completed
        with open(status_file, "w") as f:
            f.write("completed\n")
        
        # Log the output files
        log_file = job_dir / "output_files.txt"
        with open(log_file, "w") as f:
            for offset, path in sorted(output_files.items()):
                f.write(f"{offset}: {path}\n")
    
    except Exception as e:
        # Update status to failed
        with open(status_file, "w") as f:
            f.write("failed\n")
        
        # Save error information
        error_file = job_dir / "error.txt"
        with open(error_file, "w") as f:
            f.write(f"Error: {str(e)}\n\n")
            f.write("Traceback:\n")
            f.write(traceback.format_exc())


# ============================================================================
# API Endpoints
# ============================================================================

@app.get("/")
async def root():
    """Root endpoint - API information."""
    return {
        "name": "RiboStructMapper API",
        "version": "1.0.0",
        "status": "operational",
        "endpoints": {
            "submit_job": "POST /submit_job",
            "job_status": "GET /status/{job_id} (coming soon)",
            "download": "GET /download/{job_id} (coming soon)"
        }
    }


@app.post("/submit_job", response_model=JobResponse)
async def submit_job(
    background_tasks: BackgroundTasks,
    pdb_file: UploadFile = File(..., description="PDB structure file"),
    fasta_file: UploadFile = File(..., description="Genomic FASTA file"),
    gtf_file: UploadFile = File(..., description="GTF/GFF annotation file"),
    density_file: UploadFile = File(..., description="bedGraph density file"),
    offsets: str = Form(..., description="Comma-separated offset values (e.g., '0,-10,-15')")
) -> JobResponse:
    """
    Submit a new job for processing.
    
    This endpoint accepts all required input files and creates a new job.
    Files are saved to a unique job directory for later processing.
    
    Args:
        pdb_file: PDB structure file
        fasta_file: Genomic sequence in FASTA format
        gtf_file: Gene annotation in GTF/GFF format
        density_file: Ribosome density in bedGraph format
        offsets: Comma-separated offset values (e.g., "0,-10,-15")
    
    Returns:
        JobResponse with job_id and status
    
    Example:
        ```
        curl -X POST "http://localhost:8000/submit_job" \\
          -F "pdb_file=@protein.pdb" \\
          -F "fasta_file=@genome.fasta" \\
          -F "gtf_file=@annotation.gtf" \\
          -F "density_file=@density.bedgraph" \\
          -F "offsets=0,-10,-15"
        ```
    """
    try:
        # Generate unique job ID
        job_id = create_job_id()
        
        # Create job directory
        job_dir = JOBS_DIR / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        
        # Save uploaded files
        save_upload_file(pdb_file, job_dir / "input.pdb")
        save_upload_file(fasta_file, job_dir / "input.fasta")
        save_upload_file(gtf_file, job_dir / "input.gtf")
        save_upload_file(density_file, job_dir / "input.bedgraph")
        
        # Save parameters
        params_file = job_dir / "params.txt"
        with open(params_file, "w") as f:
            f.write(f"offsets={offsets}\n")
        
        # Create a status file
        status_file = job_dir / "status.txt"
        with open(status_file, "w") as f:
            f.write("uploaded\n")
        
        # Trigger background processing
        background_tasks.add_task(process_job, job_id)
        
        return JobResponse(
            job_id=job_id,
            status="queued",
            message=f"Job created and queued for processing. Job ID: {job_id}"
        )
    
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error creating job: {str(e)}"
        )


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "jobs_directory": str(JOBS_DIR),
        "jobs_directory_exists": JOBS_DIR.exists()
    }


# ============================================================================
# For development/testing
# ============================================================================

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

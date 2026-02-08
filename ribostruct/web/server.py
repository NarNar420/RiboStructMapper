"""
RiboStructMapper FastAPI Server

This module provides the web API for the RiboStructMapper pipeline.
It handles file uploads, job management, and serves processed results.
"""

import os
import uuid
import shutil
import traceback
import zipfile
import asyncio
import time
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, UploadFile, File, Form, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# Import the pipeline function
from ribostruct.cli.main import run_pipeline


# Initialize FastAPI application
app = FastAPI(
    title="RiboStructMapper API",
    description="Web API for mapping ribosome density onto PDB structures",
    version="1.0.0"
)

# Configuration
JOBS_DIR = Path("jobs")
JOBS_DIR.mkdir(exist_ok=True)

STATIC_DIR = Path(__file__).parent / "static"
STATIC_DIR.mkdir(exist_ok=True)


# ============================================================================
# Models
# ============================================================================

class JobResponse(BaseModel):
    """Response model for job submission."""
    job_id: str
    status: str
    message: Optional[str] = None


class StatusResponse(BaseModel):
    """Response model for job status."""
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
        
        # GTF is optional - check if it exists
        gtf_file_path = job_dir / "input.gtf"
        gtf_path = str(gtf_file_path) if gtf_file_path.exists() else None
        
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


def cleanup_old_jobs(max_age_hours: int = 24) -> int:
    """
    Delete job folders older than the specified age.
    
    This function prevents disk overflow by automatically cleaning up
    old completed/failed jobs.
    
    Args:
        max_age_hours: Maximum age of jobs in hours (default: 24)
    
    Returns:
        Number of jobs deleted
    """
    if not JOBS_DIR.exists():
        return 0
    
    current_time = time.time()
    max_age_seconds = max_age_hours * 3600
    deleted_count = 0
    
    for job_dir in JOBS_DIR.iterdir():
        if not job_dir.is_dir():
            continue
        
        # Get the creation/modification time of the directory
        dir_mtime = job_dir.stat().st_mtime
        age_seconds = current_time - dir_mtime
        
        if age_seconds > max_age_seconds:
            try:
                shutil.rmtree(job_dir)
                job_id = job_dir.name
                print(f"[CLEANUP] Deleted old job: {job_id} (age: {age_seconds/3600:.1f} hours)")
                deleted_count += 1
            except Exception as e:
                print(f"[CLEANUP] Error deleting job {job_dir.name}: {e}")
    
    if deleted_count > 0:
        print(f"[CLEANUP] Total jobs deleted: {deleted_count}")
    
    return deleted_count


async def background_cleanup_task():
    """
    Background task that runs cleanup periodically.
    
    This task runs every hour and cleans up jobs older than 24 hours.
    """
    while True:
        try:
            await asyncio.sleep(3600)  # Run every hour
            print("[CLEANUP] Running scheduled cleanup...")
            deleted = cleanup_old_jobs(max_age_hours=24)
            if deleted == 0:
                print("[CLEANUP] No old jobs to clean up")
        except Exception as e:
            print(f"[CLEANUP] Error in background cleanup: {e}")


# ============================================================================
# API Endpoints
# ============================================================================

@app.on_event("startup")
async def startup_event():
    """Run when the server starts."""
    print("[SERVER] Starting RiboStructMapper server...")
    print(f"[SERVER] Jobs directory: {JOBS_DIR}")
    print(f"[SERVER] Static directory: {STATIC_DIR}")
    
    # Start the background cleanup task
    asyncio.create_task(background_cleanup_task())
    print("[SERVER] Background cleanup task started (runs every hour)")


@app.get("/api")
async def root():
    """Root endpoint - API information."""
    return {
        "name": "RiboStructMapper API",
        "version": "1.0.0",
        "status": "operational",
        "endpoints": {
            "submit_job": "POST /submit_job",
            "job_status": "GET /status/{job_id}",
            "download": "GET /download/{job_id}"
        }
    }


@app.post("/submit_job", response_model=JobResponse)
async def submit_job(
    background_tasks: BackgroundTasks,
    pdb_file: UploadFile = File(..., description="PDB structure file"),
    fasta_file: UploadFile = File(..., description="FASTA file (genomic or CDS)"),
    gtf_file: Optional[UploadFile] = File(None, description="GTF/GFF annotation file (optional)"),
    density_file: UploadFile = File(..., description="bedGraph density file"),
    offsets: str = Form(..., description="Comma-separated offset values (e.g., '0,-10,-15')")
) -> JobResponse:
    """
    Submit a new job for processing.
    
    This endpoint accepts required input files and creates a new job.
    GTF file is optional - if not provided, FASTA is assumed to be CDS.
    
    Args:
        pdb_file: PDB structure file
        fasta_file: FASTA file (genomic with GTF, or CDS without GTF)
        gtf_file: Optional gene annotation in GTF/GFF format
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
        
        # GTF is optional
        if gtf_file:
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


@app.get("/status/{job_id}", response_model=StatusResponse)
async def get_job_status(job_id: str) -> StatusResponse:
    """
    Get the status of a job.
    
    Args:
        job_id: The unique job identifier
    
    Returns:
        StatusResponse with current job status
    
    Possible status values:
        - "not_found": Job does not exist
        - "uploaded": Files uploaded, awaiting processing
        - "queued": Job queued for background processing
        - "processing": Currently being processed
        - "completed": Processing finished successfully
        - "failed": Processing encountered an error
    
    Example:
        GET /status/a3b8d1b6-0b3b-4b1a-9c1a-1a2b3c4d5e6f
    """
    job_dir = JOBS_DIR / job_id
    
    # Check if job exists
    if not job_dir.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Job {job_id} not found"
        )
    
    # Read status file
    status_file = job_dir / "status.txt"
    
    if not status_file.exists():
        return StatusResponse(
            job_id=job_id,
            status="unknown",
            message="Status file not found"
        )
    
    status = status_file.read_text().strip()
    
    # Check for error file if status is failed
    message = None
    if status == "failed":
        error_file = job_dir / "error.txt"
        if error_file.exists():
            error_content = error_file.read_text()
            # Extract first line of error
            first_line = error_content.split("\n")[0]
            message = first_line
    
    return StatusResponse(
        job_id=job_id,
        status=status,
        message=message
    )


@app.get("/download/{job_id}")
async def download_results(job_id: str):
    """
    Download the results of a completed job as a ZIP file.
    
    Args:
        job_id: The unique job identifier
    
    Returns:
        FileResponse containing a ZIP archive of all output PDB files
    
    Raises:
        HTTPException 404: Job not found
        HTTPException 400: Job not completed yet
    
    Example:
        GET /download/a3b8d1b6-0b3b-4b1a-9c1a-1a2b3c4d5e6f
    """
    job_dir = JOBS_DIR / job_id
    
    # Check if job exists
    if not job_dir.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Job {job_id} not found"
        )
    
    # Check status
    status_file = job_dir / "status.txt"
    
    if not status_file.exists():
        raise HTTPException(
            status_code=400,
            detail="Job status unknown"
        )
    
    status = status_file.read_text().strip()
    
    if status != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job not completed yet. Current status: {status}"
        )
    
    # Find all output PDB files
    output_files = list(job_dir.glob("output_offset_*.pdb"))
    
    if not output_files:
        raise HTTPException(
            status_code=404,
            detail="No output files found for this job"
        )
    
    # Create ZIP file
    zip_filename = f"{job_id}_results.zip"
    zip_path = job_dir / zip_filename
    
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for output_file in output_files:
            # Add file to zip with just the filename (not full path)
            zipf.write(output_file, arcname=output_file.name)
    
    # Return the ZIP file
    return FileResponse(
        path=str(zip_path),
        media_type="application/zip",
        filename=zip_filename
    )


# ============================================================================
# Static Files - Must be last
# ============================================================================

# Mount static files (HTML frontend)
app.mount("/", StaticFiles(directory=str(STATIC_DIR), html=True), name="static")


# ============================================================================
# For development/testing
# ============================================================================

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

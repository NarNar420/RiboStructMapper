#!/usr/bin/env python3
"""
Manual test script to verify the server's background processing capability.

This script submits a job to the running server and polls the status until completion.
"""

import time
import requests
from pathlib import Path


def test_server_background_processing():
    """Submit a job and monitor its progress."""
    
    base_url = "http://127.0.0.1:8000"
    
    print("=" * 70)
    print("TESTING SERVER BACKGROUND PROCESSING")
    print("=" * 70)
    
    # Check if server is running
    try:
        response = requests.get(f"{base_url}/health")
        print(f"✓ Server is running")
        print(f"  Health: {response.json()}")
    except requests.exceptions.ConnectionError:
        print("✗ Server not running!")
        print("  Please start the server first:")
        print("  uvicorn server:app --reload")
        return
    
    # Prepare files from mock_data
    mock_dir = Path("mock_data")
    
    files = {
        "pdb_file": ("mock.pdb", open(mock_dir / "mock.pdb", "rb")),
        "fasta_file": ("mock.fasta", open(mock_dir / "mock.fasta", "rb")),
        "gtf_file": ("mock.gtf", open(mock_dir / "mock.gtf", "rb")),
        "density_file": ("mock.bedgraph", open(mock_dir / "mock.bedgraph", "rb")),
    }
    
    data = {
        "offsets": "0,-12"
    }
    
    print("\n[1] Submitting job...")
    response = requests.post(f"{base_url}/submit_job", files=files, data=data)
    
    # Close files
    for _, f in files.items():
        f[1].close()
    
    if response.status_code != 200:
        print(f"✗ Job submission failed: {response.status_code}")
        print(f"  Response: {response.text}")
        return
    
    result = response.json()
    job_id = result["job_id"]
    
    print(f"✓ Job submitted successfully")
    print(f"  Job ID: {job_id}")
    print(f"  Status: {result['status']}")
    
    # Monitor the job
    job_dir = Path("jobs") / job_id
    status_file = job_dir / "status.txt"
    
    print(f"\n[2] Monitoring job status...")
    
    max_wait = 30  # seconds
    start_time = time.time()
    
    while time.time() - start_time < max_wait:
        if status_file.exists():
            status = status_file.read_text().strip()
            elapsed = time.time() - start_time
            print(f"  [{elapsed:.1f}s] Status: {status}")
            
            if status == "completed":
                print(f"\n✓ Job completed successfully!")
                
                # List output files
                output_files = list(job_dir.glob("output_offset_*.pdb"))
                print(f"\n[3] Output files generated:")
                for f in sorted(output_files):
                    size = f.stat().st_size
                    print(f"  - {f.name} ({size} bytes)")
                
                # Check output_files.txt
                output_log = job_dir / "output_files.txt"
                if output_log.exists():
                    print(f"\n[4] Output log:")
                    print(f"  {output_log.read_text()}")
                
                break
            
            elif status == "failed":
                print(f"\n✗ Job failed!")
                error_file = job_dir / "error.txt"
                if error_file.exists():
                    print(f"\nError details:")
                    print(error_file.read_text())
                break
        
        time.sleep(1)
    else:
        print(f"\n⚠ Timeout after {max_wait}s")
    
    print("\n" + "=" * 70)


if __name__ == "__main__":
    test_server_background_processing()

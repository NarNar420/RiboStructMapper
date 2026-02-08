#!/usr/bin/env python3
"""
Complete end-to-end test for RiboStructMapper server.

Tests all endpoints: submit_job, status, and download.
Requires the server to be running.
"""

import time
import requests
import zipfile
from pathlib import Path


def main():
    """Run complete end-to-end test."""
    
    BASE_URL = "http://127.0.0.1:8000"
    
    print("="*70)
    print("RIBOSTRUCTMAPPER SERVER - COMPLETE END-TO-END TEST")
    print("="*70)
    
    # Step 1: Check server is running
    print("\n[1/5] Checking server health...")
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=2)
        response.raise_for_status()
        health = response.json()
        print(f"✓ Server is healthy")
        print(f"  Jobs directory: {health['jobs_directory']}")
    except requests.exceptions.ConnectionError:
        print("✗ Server is not running!")
        print(f"\nPlease start the server first:")
        print(f"  uvicorn server:app --reload")
        print(f"\nThen run this test again.")
        return False
    except Exception as e:
        print(f"✗ Health check failed: {e}")
        return False
    
    # Step 2: Submit a job
    print(f"\n[2/5] Submitting a test job...")
    
    mock_dir = Path("mock_data")
    if not mock_dir.exists():
        print(f"✗ mock_data directory not found")
        return False
    
    try:
        with open(mock_dir / "mock.pdb", "rb") as pdb, \
             open(mock_dir / "mock.fasta", "rb") as fasta, \
             open(mock_dir / "mock.gtf", "rb") as gtf, \
             open(mock_dir / "mock.bedgraph", "rb") as bedgraph:
            
            files = {
                "pdb_file": ("mock.pdb", pdb),
                "fasta_file": ("mock.fasta", fasta),
                "gtf_file": ("mock.gtf", gtf),
                "density_file": ("mock.bedgraph", bedgraph),
            }
            
            data = {"offsets": "0,-12"}
            
            response = requests.post(f"{BASE_URL}/submit_job", files=files, data=data)
            response.raise_for_status()
            
        result = response.json()
        job_id = result["job_id"]
        
        print(f"✓ Job submitted successfully")
        print(f"  Job ID: {job_id}")
        print(f"  Status: {result['status']}")
        
    except Exception as e:
        print(f"✗ Job submission failed: {e}")
        return False
    
    # Step 3: Poll status until completed
    print(f"\n[3/5] Monitoring job status...")
    
    max_wait = 30  # seconds
    start_time = time.time()
    last_status = None
    
    while time.time() - start_time < max_wait:
        try:
            response = requests.get(f"{BASE_URL}/status/{job_id}")
            response.raise_for_status()
            
            status_data = response.json()
            current_status = status_data["status"]
            
            # Print status changes
            if current_status != last_status:
                elapsed = time.time() - start_time
                print(f"  [{elapsed:.1f}s] Status: {current_status}")
                last_status = current_status
            
            if current_status == "completed":
                print(f"✓ Job completed successfully")
                break
            elif current_status == "failed":
                error_msg = status_data.get("message", "Unknown error")
                print(f"✗ Job failed: {error_msg}")
                return False
            
            time.sleep(0.5)
            
        except Exception as e:
            print(f"✗ Status check failed: {e}")
            return False
    else:
        print(f"✗ Timeout waiting for job completion ({max_wait}s)")
        return False
    
    # Step 4: Test status endpoint for non-existent job
    print(f"\n[4/5] Testing error handling...")
    
    try:
        response = requests.get(f"{BASE_URL}/status/nonexistent-id")
        
        if response.status_code == 404:
            print(f"✓ Status endpoint returns 404 for non-existent job")
        else:
            print(f"✗ Expected 404, got {response.status_code}")
            return False
            
    except Exception as e:
        print(f"✗ Error handling test failed: {e}")
        return False
    
    # Step 5: Download results
    print(f"\n[5/5] Downloading results...")
    
    try:
        response = requests.get(f"{BASE_URL}/download/{job_id}")
        response.raise_for_status()
        
        if response.headers.get("content-type") != "application/zip":
            print(f"✗ Expected ZIP file, got {response.headers.get('content-type')}")
            return False
        
        # Save the ZIP file
        zip_path = Path(f"test_download_{job_id[:8]}.zip")
        with open(zip_path, "wb") as f:
            f.write(response.content)
        
        file_size = zip_path.stat().st_size
        print(f"✓ Downloaded ZIP file: {zip_path} ({file_size} bytes)")
        
        # Verify ZIP contents
        with zipfile.ZipFile(zip_path, 'r') as zipf:
            files_in_zip = zipf.namelist()
            print(f"✓ ZIP contains {len(files_in_zip)} files:")
            for filename in sorted(files_in_zip):
                print(f"    - {filename}")
        
        # Clean up
        zip_path.unlink()
        print(f"✓ Test cleanup complete")
        
    except Exception as e:
        print(f"✗ Download failed: {e}")
        return False
    
    # Success!
    print(f"\n" + "="*70)
    print("✅ ALL TESTS PASSED!")
    print("="*70)
    print(f"\nThe server is working correctly with all endpoints:")
    print(f"  • POST /submit_job - ✓")
    print(f"  • GET /status/{{job_id}} - ✓")
    print(f"  • GET /download/{{job_id}} - ✓")
    print(f"\nYou can now access the server at:")
    print(f"  API Docs: {BASE_URL}/docs")
    print(f"  Health: {BASE_URL}/health")
    print(f"="*70 + "\n")
    
    return True


if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)

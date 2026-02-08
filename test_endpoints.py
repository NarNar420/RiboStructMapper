#!/usr/bin/env python3
"""
Comprehensive test for status and download endpoints.

This script tests the complete workflow including:
1. Job submission
2. Status polling
3. Result download
"""

import time
import requests
import zipfile
from pathlib import Path


def test_status_and_download():
    """Test the complete workflow with status and download endpoints."""
    
    base_url = "http://127.0.0.1:8000"
    
    print("=" * 70)
    print("TESTING STATUS AND DOWNLOAD ENDPOINTS")
    print("=" * 70)
    
    # Check if server is running
    try:
        response = requests.get(f"{base_url}/health")
        print(f"✓ Server is running")
    except requests.exceptions.ConnectionError:
        print("✗ Server not running!")
        print("  Start with: uvicorn server:app --reload")
        return False
    
    # Test 1: Status endpoint for non-existent job
    print("\n[1] Testing status endpoint for non-existent job...")
    response = requests.get(f"{base_url}/status/nonexistent-job-id")
    assert response.status_code == 404, "Should return 404 for non-existent job"
    print(f"✓ Returns 404 for non-existent job")
    
    # Test 2: Submit a job
    print("\n[2] Submitting a test job...")
    mock_dir = Path("mock_data")
    
    files = {
        "pdb_file": ("mock.pdb", open(mock_dir / "mock.pdb", "rb")),
        "fasta_file": ("mock.fasta", open(mock_dir / "mock.fasta", "rb")),
        "gtf_file": ("mock.gtf", open(mock_dir / "mock.gtf", "rb")),
        "density_file": ("mock.bedgraph", open(mock_dir / "mock.bedgraph", "rb")),
    }
    
    data = {"offsets": "0,-12"}
    
    response = requests.post(f"{base_url}/submit_job", files=files, data=data)
    
    # Close files
    for _, f in files.items():
        f[1].close()
    
    assert response.status_code == 200
    result = response.json()
    job_id = result["job_id"]
    
    print(f"✓ Job submitted: {job_id}")
    print(f"  Initial status: {result['status']}")
    
    # Test 3: Poll status until completed
    print(f"\n[3] Polling job status...")
    
    max_wait = 30
    start_time = time.time()
    status_history = []
    
    while time.time() - start_time < max_wait:
        response = requests.get(f"{base_url}/status/{job_id}")
        assert response.status_code == 200
        
        status_data = response.json()
        current_status = status_data["status"]
        
        # Track status changes
        if not status_history or status_history[-1] != current_status:
            elapsed = time.time() - start_time
            print(f"  [{elapsed:.1f}s] Status: {current_status}")
            status_history.append(current_status)
        
        if current_status == "completed":
            print(f"✓ Job completed successfully")
            break
        elif current_status == "failed":
            print(f"✗ Job failed")
            print(f"  Error: {status_data.get('message', 'Unknown error')}")
            return False
        
        time.sleep(0.5)
    else:
        print(f"✗ Timeout waiting for job completion")
        return False
    
    # Test 4: Try to download before completion would fail (already completed in this test)
    # So we test the download endpoint directly
    print(f"\n[4] Testing download endpoint...")
    
    response = requests.get(f"{base_url}/download/{job_id}")
    assert response.status_code == 200, f"Download should succeed, got {response.status_code}"
    assert response.headers["content-type"] == "application/zip"
    
    print(f"✓ Download endpoint returns ZIP file")
    
    # Save and verify ZIP contents
    zip_path = Path(f"test_download_{job_id}.zip")
    with open(zip_path, "wb") as f:
        f.write(response.content)
    
    print(f"✓ ZIP file saved: {zip_path} ({zip_path.stat().st_size} bytes)")
    
    # Extract and verify contents
    with zipfile.ZipFile(zip_path, 'r') as zipf:
        file_list = zipf.namelist()
        print(f"✓ ZIP contains {len(file_list)} files:")
        for filename in file_list:
            print(f"    - {filename}")
    
    # Check that we have the expected output files
    expected_files = {"output_offset_0.pdb", "output_offset_-12.pdb"}
    actual_files = set(file_list)
    
    assert expected_files == actual_files, f"Expected {expected_files}, got {actual_files}"
    print(f"✓ All expected output files present in ZIP")
    
    # Test 5: Download non-existent job
    print(f"\n[5] Testing download for non-existent job...")
    response = requests.get(f"{base_url}/download/nonexistent-job-id")
    assert response.status_code == 404
    print(f"✓ Returns 404 for non-existent job")
    
    # Cleanup
    print(f"\n[6] Cleaning up test files...")
    if zip_path.exists():
        zip_path.unlink()
    print(f"✓ Cleanup complete")
    
    print("\n" + "=" * 70)
    print("ALL TESTS PASSED ✓")
    print("=" * 70)
    
    return True


if __name__ == "__main__":
    import sys
    success = test_status_and_download()
    sys.exit(0 if success else 1)

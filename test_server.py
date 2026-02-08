#!/usr/bin/env python3
"""
Test suite for the FastAPI server.

This script tests the server's file upload and job management functionality
using FastAPI's TestClient.
"""

import os
import sys
import shutil
from pathlib import Path

from fastapi.testclient import TestClient
from server import app, JOBS_DIR


# Initialize test client
client = TestClient(app)


def cleanup_jobs_dir():
    """Clean up the jobs directory before/after tests."""
    if JOBS_DIR.exists():
        shutil.rmtree(JOBS_DIR)
    JOBS_DIR.mkdir(exist_ok=True)


def test_root_endpoint():
    """Test the root endpoint returns API information."""
    print("=" * 60)
    print("Test 1: Root endpoint")
    print("=" * 60)
    
    try:
        response = client.get("/")
        
        assert response.status_code == 200, f"Expected 200, got {response.status_code}"
        data = response.json()
        
        assert "name" in data, "Response should contain 'name'"
        assert "RiboStructMapper" in data["name"], "Name should contain 'RiboStructMapper'"
        
        print(f"✓ Root endpoint working")
        print(f"  API Name: {data['name']}")
        print(f"  Version: {data['version']}")
        print("\n✓✓✓ Test 1 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 1 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_health_check():
    """Test the health check endpoint."""
    print("=" * 60)
    print("Test 2: Health check endpoint")
    print("=" * 60)
    
    try:
        response = client.get("/health")
        
        assert response.status_code == 200
        data = response.json()
        
        assert data["status"] == "healthy"
        assert data["jobs_directory_exists"] == True
        
        print(f"✓ Health check passed")
        print(f"  Status: {data['status']}")
        print(f"  Jobs directory: {data['jobs_directory']}")
        print("\n✓✓✓ Test 2 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 2 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_submit_job():
    """Test job submission with file uploads."""
    print("=" * 60)
    print("Test 3: Submit job with file uploads")
    print("=" * 60)
    
    cleanup_jobs_dir()
    
    try:
        # Prepare test files from mock_data
        mock_dir = Path("mock_data")
        
        files = {
            "pdb_file": ("mock.pdb", open(mock_dir / "mock.pdb", "rb"), "chemical/x-pdb"),
            "fasta_file": ("mock.fasta", open(mock_dir / "mock.fasta", "rb"), "text/plain"),
            "gtf_file": ("mock.gtf", open(mock_dir / "mock.gtf", "rb"), "text/plain"),
            "density_file": ("mock.bedgraph", open(mock_dir / "mock.bedgraph", "rb"), "text/plain"),
        }
        
        data = {
            "offsets": "0,-10,-15"
        }
        
        # Submit the job
        response = client.post("/submit_job", files=files, data=data)
        
        # Close all file handles
        for _, f in files.items():
            f[1].close()
        
        print(f"  Response status: {response.status_code}")
        
        assert response.status_code == 200, f"Expected 200, got {response.status_code}"
        
        response_data = response.json()
        print(f"  Response: {response_data}")
        
        assert "job_id" in response_data, "Response should contain 'job_id'"
        assert "status" in response_data, "Response should contain 'status'"
        assert response_data["status"] == "uploaded", f"Status should be 'uploaded', got {response_data['status']}"
        
        job_id = response_data["job_id"]
        print(f"✓ Job created with ID: {job_id}")
        
        # Verify job directory was created
        job_dir = JOBS_DIR / job_id
        assert job_dir.exists(), f"Job directory should exist: {job_dir}"
        print(f"✓ Job directory created: {job_dir}")
        
        # Verify all files were saved
        expected_files = ["input.pdb", "input.fasta", "input.gtf", "input.bedgraph", "params.txt", "status.txt"]
        for filename in expected_files:
            file_path = job_dir / filename
            assert file_path.exists(), f"File should exist: {file_path}"
        
        print(f"✓ All expected files saved:")
        for filename in expected_files:
            file_path = job_dir / filename
            size = file_path.stat().st_size
            print(f"    {filename}: {size} bytes")
        
        # Verify params.txt content
        params_content = (job_dir / "params.txt").read_text()
        assert "offsets=0,-10,-15" in params_content, "params.txt should contain offsets"
        print(f"✓ Parameters saved correctly")
        
        # Verify status.txt content
        status_content = (job_dir / "status.txt").read_text().strip()
        assert status_content == "uploaded", "status.txt should contain 'uploaded'"
        print(f"✓ Status file created correctly")
        
        print("\n✓✓✓ Test 3 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 3 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_multiple_job_submissions():
    """Test that multiple job submissions create separate directories."""
    print("=" * 60)
    print("Test 4: Multiple job submissions")
    print("=" * 60)
    
    cleanup_jobs_dir()
    
    try:
        job_ids = []
        
        # Submit 3 jobs
        for i in range(3):
            mock_dir = Path("mock_data")
            
            files = {
                "pdb_file": ("mock.pdb", open(mock_dir / "mock.pdb", "rb")),
                "fasta_file": ("mock.fasta", open(mock_dir / "mock.fasta", "rb")),
                "gtf_file": ("mock.gtf", open(mock_dir / "mock.gtf", "rb")),
                "density_file": ("mock.bedgraph", open(mock_dir / "mock.bedgraph", "rb")),
            }
            
            data = {"offsets": f"0,-{i*10}"}
            
            response = client.post("/submit_job", files=files, data=data)
            
            # Close files
            for _, f in files.items():
                f[1].close()
            
            assert response.status_code == 200
            job_id = response.json()["job_id"]
            job_ids.append(job_id)
        
        print(f"✓ Created {len(job_ids)} jobs")
        
        # Verify all job IDs are unique
        assert len(job_ids) == len(set(job_ids)), "All job IDs should be unique"
        print(f"✓ All job IDs are unique")
        
        # Verify all directories exist
        for job_id in job_ids:
            job_dir = JOBS_DIR / job_id
            assert job_dir.exists(), f"Job directory should exist: {job_dir}"
        
        print(f"✓ All job directories created")
        for i, job_id in enumerate(job_ids):
            print(f"    Job {i+1}: {job_id}")
        
        print("\n✓✓✓ Test 4 PASSED ✓✓✓\n")
        return True
        
    except Exception as e:
        print(f"✗ Test 4 FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("SERVER TEST SUITE")
    print("=" * 60 + "\n")
    
    # Make sure we're in the right directory
    if not Path("mock_data").exists():
        print("Error: mock_data directory not found. Run from project root.")
        sys.exit(1)
    
    results = []
    results.append(("Test 1: Root endpoint", test_root_endpoint()))
    results.append(("Test 2: Health check", test_health_check()))
    results.append(("Test 3: Submit job", test_submit_job()))
    results.append(("Test 4: Multiple jobs", test_multiple_job_submissions()))
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    print("=" * 60 + "\n")
    
    # Cleanup
    print("Cleaning up test jobs directory...")
    cleanup_jobs_dir()
    print("✓ Cleanup complete\n")
    
    return all(result for _, result in results)


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

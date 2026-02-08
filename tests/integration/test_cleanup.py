#!/usr/bin/env python3
"""
Test script for automated job cleanup functionality.

This script tests that old jobs are correctly identified and deleted.
"""

import os
import time
import shutil
from pathlib import Path

# Import the cleanup function from server
from ribostruct.web.server import cleanup_old_jobs, JOBS_DIR



def test_cleanup_old_jobs():
    """Test that jobs older than 24 hours are deleted."""
    print("=" * 60)
    print("TEST: Cleanup Old Jobs")
    print("=" * 60)
    
    # Ensure jobs directory exists
    JOBS_DIR.mkdir(exist_ok=True)
    
    # Create test job directories
    test_jobs = []
    
    # Job 1: Recent job (should NOT be deleted)
    recent_job_id = "test-job-recent"
    recent_job_dir = JOBS_DIR / recent_job_id
    recent_job_dir.mkdir(exist_ok=True)
    (recent_job_dir / "test.txt").write_text("recent job")
    test_jobs.append((recent_job_id, recent_job_dir, False))  # False = should not delete
    
    # Job 2: Old job - 25 hours ago (should be deleted)
    old_job_id = "test-job-old-25h"
    old_job_dir = JOBS_DIR / old_job_id
    old_job_dir.mkdir(exist_ok=True)
    (old_job_dir / "test.txt").write_text("old job 25h")
    
    # Set modification time to 25 hours ago
    hours_ago_25 = time.time() - (25 * 3600)
    os.utime(old_job_dir, (hours_ago_25, hours_ago_25))
    test_jobs.append((old_job_id, old_job_dir, True))  # True = should delete
    
    # Job 3: Very old job - 48 hours ago (should be deleted)
    very_old_job_id = "test-job-old-48h"
    very_old_job_dir = JOBS_DIR / very_old_job_id
    very_old_job_dir.mkdir(exist_ok=True)
    (very_old_job_dir / "test.txt").write_text("old job 48h")
    
    # Set modification time to 48 hours ago
    hours_ago_48 = time.time() - (48 * 3600)
    os.utime(very_old_job_dir, (hours_ago_48, hours_ago_48))
    test_jobs.append((very_old_job_id, very_old_job_dir, True))  # True = should delete
    
    # Job 4: Borderline job - 23 hours ago (should NOT be deleted)
    borderline_job_id = "test-job-borderline-23h"
    borderline_job_dir = JOBS_DIR / borderline_job_id
    borderline_job_dir.mkdir(exist_ok=True)
    (borderline_job_dir / "test.txt").write_text("borderline job")
    
    # Set modification time to 23 hours ago
    hours_ago_23 = time.time() - (23 * 3600)
    os.utime(borderline_job_dir, (hours_ago_23, hours_ago_23))
    test_jobs.append((borderline_job_id, borderline_job_dir, False))  # False = should not delete
    
    print(f"\nCreated {len(test_jobs)} test job directories:")
    for job_id, job_dir, should_delete in test_jobs:
        age_hours = (time.time() - job_dir.stat().st_mtime) / 3600
        print(f"  - {job_id}: {age_hours:.1f} hours old (should {'DELETE' if should_delete else 'KEEP'})")
    
    # Run cleanup
    print(f"\nRunning cleanup (max age: 24 hours)...")
    deleted_count = cleanup_old_jobs(max_age_hours=24)
    
    print(f"\nVerifying results...")
    
    # Verify each job
    all_passed = True
    for job_id, job_dir, should_delete in test_jobs:
        exists = job_dir.exists()
        
        if should_delete:
            if exists:
                print(f"  ✗ FAIL: {job_id} should have been deleted but still exists")
                all_passed = False
            else:
                print(f"  ✓ PASS: {job_id} was correctly deleted")
        else:
            if exists:
                print(f"  ✓ PASS: {job_id} was correctly kept")
            else:
                print(f"  ✗ FAIL: {job_id} should have been kept but was deleted")
                all_passed = False
    
    # Verify count
    expected_deletions = sum(1 for _, _, should_delete in test_jobs if should_delete)
    if deleted_count == expected_deletions:
        print(f"\n✓ Deletion count correct: {deleted_count}/{expected_deletions}")
    else:
        print(f"\n✗ Deletion count mismatch: {deleted_count} (expected {expected_deletions})")
        all_passed = False
    
    # Cleanup remaining test jobs
    print(f"\nCleaning up remaining test jobs...")
    for job_id, job_dir, _ in test_jobs:
        if job_dir.exists():
            shutil.rmtree(job_dir)
    
    print("\n" + "=" * 60)
    if all_passed:
        print("✅ ALL CLEANUP TESTS PASSED")
    else:
        print("❌ SOME CLEANUP TESTS FAILED")
    print("=" * 60 + "\n")
    
    return all_passed


def test_cleanup_with_custom_age():
    """Test cleanup with custom age threshold."""
    print("=" * 60)
    print("TEST: Cleanup with Custom Age (12 hours)")
    print("=" * 60)
    
    # Create a job 15 hours old
    test_job_id = "test-job-15h"
    test_job_dir = JOBS_DIR / test_job_id
    test_job_dir.mkdir(exist_ok=True)
    (test_job_dir / "test.txt").write_text("test job")
    
    hours_ago_15 = time.time() - (15 * 3600)
    os.utime(test_job_dir, (hours_ago_15, hours_ago_15))
    
    age_hours = (time.time() - test_job_dir.stat().st_mtime) / 3600
    print(f"\nCreated test job: {test_job_id} ({age_hours:.1f} hours old)")
    
    # Run cleanup with 12-hour threshold
    print(f"Running cleanup with max_age_hours=12...")
    deleted_count = cleanup_old_jobs(max_age_hours=12)
    
    if not test_job_dir.exists() and deleted_count == 1:
        print(f"✓ Job older than 12 hours was correctly deleted")
        result = True
    else:
        print(f"✗ Job deletion failed")
        result = False
        if test_job_dir.exists():
            shutil.rmtree(test_job_dir)
    
    print("\n" + "=" * 60)
    if result:
        print("✅ CUSTOM AGE TEST PASSED")
    else:
        print("❌ CUSTOM AGE TEST FAILED")
    print("=" * 60 + "\n")
    
    return result


def test_cleanup_empty_directory():
    """Test cleanup when jobs directory is empty."""
    print("=" * 60)
    print("TEST: Cleanup Empty Directory")
    print("=" * 60)
    
    # Ensure directory is empty
    if JOBS_DIR.exists():
        for item in JOBS_DIR.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
    
    print("\nRunning cleanup on empty directory...")
    deleted_count = cleanup_old_jobs(max_age_hours=24)
    
    if deleted_count == 0:
        print("✓ Cleanup correctly returned 0 for empty directory")
        result = True
    else:
        print(f"✗ Cleanup returned {deleted_count} (expected 0)")
        result = False
    
    print("\n" + "=" * 60)
    if result:
        print("✅ EMPTY DIRECTORY TEST PASSED")
    else:
        print("❌ EMPTY DIRECTORY TEST FAILED")
    print("=" * 60 + "\n")
    
    return result


def main():
    """Run all cleanup tests."""
    print("\n" + "=" * 60)
    print("AUTOMATED JOB CLEANUP TEST SUITE")
    print("=" * 60 + "\n")
    
    results = []
    results.append(("Old Jobs Cleanup", test_cleanup_old_jobs()))
    results.append(("Custom Age Threshold", test_cleanup_with_custom_age()))
    results.append(("Empty Directory", test_cleanup_empty_directory()))
    
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
    
    return all(result for _, result in results)


if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)

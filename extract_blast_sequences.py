#!/usr/bin/env python3
"""
Script to extract sequences from BLAST database for eDNA pipeline processing.
"""

import subprocess
import sys
import os
from pathlib import Path
import logging

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def extract_sequences_from_blast_db(db_path: str, db_name: str, output_path: str) -> bool:
    """
    Extract sequences from BLAST database using blastdbcmd.
    
    Args:
        db_path: Path to BLAST database directory
        db_name: Name of the database (without extension)
        output_path: Output FASTA file path
        
    Returns:
        True if successful, False otherwise
    """
    logger = logging.getLogger(__name__)
    
    # Full database path
    full_db_path = os.path.join(db_path, db_name)
    
    logger.info(f"Extracting sequences from {full_db_path}")
    logger.info(f"Output file: {output_path}")
    
    try:
        # Create output directory if it doesn't exist
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        
        # Command to extract all sequences from BLAST database
        cmd = [
            "blastdbcmd",
            "-db", full_db_path,
            "-entry", "all",
            "-out", output_path,
            "-outfmt", "%f"  # FASTA format
        ]
        
        logger.info(f"Running command: {' '.join(cmd)}")
        
        # Run the command
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True,
            timeout=3600  # 1 hour timeout for large databases
        )
        
        if result.returncode == 0:
            logger.info(f"Successfully extracted sequences to {output_path}")
            
            # Check if file was created and has content
            if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                logger.info(f"Output file size: {os.path.getsize(output_path)} bytes")
                return True
            else:
                logger.error("Output file was not created or is empty")
                return False
        else:
            logger.error(f"blastdbcmd failed with return code: {result.returncode}")
            logger.error(f"Error output: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error("Command timed out after 1 hour")
        return False
    except FileNotFoundError:
        logger.error("blastdbcmd command not found. Please ensure BLAST+ is installed and in PATH")
        return False
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        return False

def count_sequences_in_fasta(fasta_path: str) -> int:
    """
    Count the number of sequences in a FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Number of sequences
    """
    count = 0
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception as e:
        logging.error(f"Error counting sequences: {e}")
        return 0
    
    return count

def main():
    """Main execution function."""
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Configuration
    blast_db_path = "C:\\blastdb"
    db_name = "16S_ribosomal_RNA"
    output_dir = "data\\raw"
    output_file = "16S_ribosomal_sequences.fasta"
    output_path = os.path.join(output_dir, output_file)
    
    logger.info("Starting BLAST database sequence extraction")
    logger.info(f"Database path: {blast_db_path}")
    logger.info(f"Database name: {db_name}")
    
    # Check if BLAST database files exist
    db_files = [f"{db_name}.nhr", f"{db_name}.nin", f"{db_name}.nsq"]
    missing_files = []
    
    for db_file in db_files:
        full_path = os.path.join(blast_db_path, db_file)
        if not os.path.exists(full_path):
            missing_files.append(db_file)
    
    if missing_files:
        logger.error(f"Missing BLAST database files: {missing_files}")
        logger.error("Please ensure the BLAST database is properly installed")
        sys.exit(1)
    
    # Extract sequences
    logger.info("Starting sequence extraction...")
    success = extract_sequences_from_blast_db(blast_db_path, db_name, output_path)
    
    if success:
        # Count sequences
        seq_count = count_sequences_in_fasta(output_path)
        logger.info(f"Extraction completed successfully!")
        logger.info(f"Total sequences extracted: {seq_count:,}")
        logger.info(f"Output file: {os.path.abspath(output_path)}")
        
        print("\n" + "="*60)
        print("       BLAST DATABASE EXTRACTION COMPLETED")
        print("="*60)
        print(f"Database: {db_name}")
        print(f"Total sequences: {seq_count:,}")
        print(f"Output file: {os.path.abspath(output_path)}")
        print(f"File size: {os.path.getsize(output_path):,} bytes")
        print("="*60)
        print("\nReady for eDNA pipeline training!")
        print("Run: python main.py --data data\\raw --blast-db C:\\blastdb")
        print("="*60)
        
    else:
        logger.error("Sequence extraction failed")
        print("\nSequence extraction failed. Please check the logs above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
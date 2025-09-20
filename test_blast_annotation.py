#!/usr/bin/env python3
"""
Test script to check BLAST annotation capabilities with real taxonomic names.
"""

import sys
import os
sys.path.append('src')

from src.utils.blast_taxonomy import BLASTTaxonomyAnnotator
from src.preprocessing.data_loader import eDNADataLoader
import pandas as pd

def test_blast_annotation():
    print("Testing BLAST annotation with sample sequences...")
    
    # Load first 50 sequences for testing
    data_loader = eDNADataLoader()
    sequences = data_loader.load_fasta_file("data/raw/16S_ribosomal_sequences.fasta")
    
    # Take first 50 sequences for testing
    sample_sequences = sequences[:50]
    print(f"Testing with {len(sample_sequences)} sequences")
    
    # Create sequence dictionary
    seq_dict = {f"seq_{i}": str(seq.seq) for i, seq in enumerate(sample_sequences)}
    
    # Initialize BLAST annotator
    try:
        blast_annotator = BLASTTaxonomyAnnotator("C:\\blastdb")
        
        # Perform BLAST annotation
        print("Running BLAST annotation...")
        blast_results = blast_annotator.annotate_sequences(seq_dict)
        
        if not blast_results.empty:
            print(f"Found {len(blast_results)} BLAST hits!")
            print("\nSample BLAST results:")
            print(blast_results[['qseqid', 'sseqid', 'pident', 'organism']].head(10))
            
            # First get consensus taxonomy to see what's available
            consensus_taxonomy = blast_annotator.get_consensus_taxonomy(blast_results, 80.0)
            print(f"\nConsensus taxonomy results ({len(consensus_taxonomy)} sequences):")
            for seq_id, taxonomy in list(consensus_taxonomy.items())[:15]:
                print(f"  {seq_id}: {taxonomy}")
            
            # Get taxonomic labels
            taxonomic_labels, detailed_results = blast_annotator.create_training_labels(
                seq_dict, min_confidence=80.0  # Lower threshold
            )
            
            print(f"\nCreated {len(taxonomic_labels)} taxonomic labels:")
            for seq_id, taxonomy in list(taxonomic_labels.items())[:10]:
                print(f"  {seq_id}: {taxonomy}")
            
            # Show unique taxonomic groups found
            unique_taxa = set(taxonomic_labels.values())
            print(f"\nUnique taxonomic groups found: {len(unique_taxa)}")
            for taxon in sorted(unique_taxa):
                count = sum(1 for t in taxonomic_labels.values() if t == taxon)
                print(f"  {taxon}: {count} sequences")
                
        else:
            print("No BLAST hits found!")
            
    except Exception as e:
        print(f"BLAST annotation failed: {e}")
        print("This might be due to BLAST not being available or database issues")

if __name__ == "__main__":
    test_blast_annotation()
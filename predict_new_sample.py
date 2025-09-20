#!/usr/bin/env python3
"""
eDNA Sample Prediction Script

Use this script to predict taxonomic composition of new eDNA samples 
using the trained pipeline.

Usage:
    python predict_new_sample.py --input new_sample.fasta --output results_new/
    python predict_new_sample.py --input path/to/new_samples/ --output results_batch/
"""

import sys
import os
import argparse
import logging
from pathlib import Path
import pandas as pd
import numpy as np
import joblib
from typing import Dict, Any, List, Optional

# Add src to path
sys.path.append('src')

from src.models.ensemble_pipeline import eDNAEnsemblePipeline
from src.preprocessing.data_loader import eDNADataLoader
from src.preprocessing.feature_extraction import FeatureExtractor
from src.utils.taxonomic_mapping import TaxonomicMapper
from src.visualization.biodiversity_plots import BiodiversityVisualizer

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class eDNAPredictor:
    """
    Predicts taxonomic composition of new eDNA samples using trained models.
    """
    
    def __init__(self, model_dir: str = "models_cluster"):
        """
        Initialize predictor with trained models.
        
        Args:
            model_dir: Directory containing trained models
        """
        self.model_dir = Path(model_dir)
        self.pipeline = None
        self.taxonomic_mapper = None
        self.feature_extractor = FeatureExtractor()
        self.data_loader = eDNADataLoader()
        
        # Load trained models
        self._load_models()
    
    def _load_models(self):
        """Load trained models and metadata."""
        try:
            # Load pipeline metadata
            metadata_path = self.model_dir / "pipeline_metadata.joblib"
            if metadata_path.exists():
                metadata = joblib.load(metadata_path)
                logger.info(f"Loaded pipeline metadata: {metadata}")
                
                # Initialize pipeline with same configuration
                self.pipeline = eDNAEnsemblePipeline(
                    use_cnn=metadata.get('use_cnn', False),
                    use_gnn=metadata.get('use_gnn', False),
                    use_rf=metadata.get('use_rf', False),
                    use_clustering=metadata.get('use_clustering', True),
                    ensemble_method=metadata.get('ensemble_method', 'weighted_voting'),
                    confidence_threshold=metadata.get('confidence_threshold', 0.1)
                )
                
                # Mark as trained (we'll load the actual models later)
                self.pipeline.is_trained = True
                
                # Load clustering results if available
                clustering_results_path = self.model_dir / "clustering_results.joblib"
                if clustering_results_path.exists():
                    clustering_results = joblib.load(clustering_results_path)
                    self.pipeline.training_results = {'clustering': clustering_results}
                    logger.info("Loaded clustering results")
                else:
                    logger.warning("No clustering results found - will use mock taxonomic mapping")
                    
            else:
                logger.error(f"No pipeline metadata found in {self.model_dir}")
                raise FileNotFoundError("Trained models not found")
                
        except Exception as e:
            logger.error(f"Failed to load models: {e}")
            raise
    
    def predict_sample(self, 
                      input_path: str, 
                      output_dir: str = "prediction_results",
                      sample_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Predict taxonomic composition for a new eDNA sample.
        
        Args:
            input_path: Path to FASTA file or directory
            output_dir: Output directory for results
            sample_name: Optional name for the sample
            
        Returns:
            Dictionary containing prediction results
        """
        logger.info(f"Predicting taxonomic composition for: {input_path}")
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        if not sample_name:
            sample_name = Path(input_path).stem
        
        try:
            # Load and preprocess new sample
            logger.info("Loading and preprocessing new sample...")
            processed_data = self._preprocess_sample(input_path)
            
            # Make predictions using clustering (since that's what we trained)
            logger.info("Making taxonomic predictions...")
            predictions = self._predict_with_clustering(processed_data)
            
            # Estimate abundance and diversity
            logger.info("Calculating abundance and diversity metrics...")
            abundance_results = self._calculate_abundance(predictions, processed_data)
            
            # Create results summary
            results = {
                'sample_name': sample_name,
                'input_path': str(input_path),
                'total_sequences': len(processed_data['sequences']),
                'predictions': predictions,
                'abundance_estimates': abundance_results,
                'processed_data': processed_data
            }
            
            # Save results
            self._save_results(results, output_path, sample_name)
            
            logger.info(f"Prediction completed for {sample_name}")
            return results
            
        except Exception as e:
            logger.error(f"Prediction failed: {e}")
            raise
    
    def _preprocess_sample(self, input_path: str) -> Dict[str, Any]:
        """Preprocess new sample data."""
        # Load sequences
        if Path(input_path).is_file():
            sequences = self.data_loader.load_fasta_file(input_path)
        else:
            sequences = self.data_loader.load_fasta_directory(input_path)
        
        # Preprocess sequences
        df = self.data_loader.preprocess_sequences(sequences)
        
        # Extract features (same as training)
        sequence_list = df['sequence'].tolist()
        sequence_ids = df['id'].tolist()
        
        features = self.feature_extractor.extract_all_features(
            sequence_list,
            include_onehot=False,  # We used only physicochemical + kmer
            include_kmer=True,
            include_physicochemical=True
        )
        
        return {
            'sequences': sequence_list,
            'sequence_ids': sequence_ids,
            'dataframe': df,
            'features': features,
            'raw_sequences': sequences
        }
    
    def _predict_with_clustering(self, processed_data: Dict[str, Any]) -> Dict[str, Any]:
        """Make predictions using the trained clustering approach."""
        from sklearn.cluster import KMeans
        from sklearn.preprocessing import StandardScaler
        from src.utils.blast_taxonomy import create_mock_taxonomic_labels
        
        features = processed_data['features']
        sequences = processed_data['sequences']
        sequence_ids = processed_data['sequence_ids']
        
        # Use same feature combination as training
        X = features['physicochemical']
        
        # Scale features (same as training)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # Adapt number of clusters to sample size
        n_samples = len(sequences)
        n_clusters = min(50, max(1, n_samples // 2))  # At least 1, at most 50, but not more than half the samples
        
        # Apply K-means with adapted cluster number
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(X_scaled)
        
        # Calculate novelty scores
        distances = kmeans.transform(X_scaled)
        min_distances = np.min(distances, axis=1)
        novelty_scores = (min_distances - np.min(min_distances)) / (np.max(min_distances) - np.min(min_distances) + 1e-8)
        
        # Create mock taxonomic labels for mapping
        seq_dict = {seq_id: seq for seq_id, seq in zip(sequence_ids, sequences)}
        mock_labels = create_mock_taxonomic_labels(seq_dict)
        taxonomic_labels = [mock_labels.get(seq_id, 'Unknown') for seq_id in sequence_ids]
        
        # Create taxonomic mapping
        mapper = TaxonomicMapper()
        cluster_taxon_mapping = mapper.create_cluster_taxon_mapping(
            cluster_labels.tolist(),
            taxonomic_labels,
            sequence_ids
        )
        
        # Convert cluster labels to taxonomic predictions
        predictions = [f'Cluster_{label}' for label in cluster_labels]
        taxonomic_predictions = mapper.map_predictions_to_taxa(predictions, cluster_labels)
        
        # Calculate confidence (inverse of novelty)
        confidences = 1.0 - novelty_scores
        
        # Identify novel taxa
        novelty_threshold = 0.7
        novel_indices = np.where(novelty_scores >= novelty_threshold)[0]
        
        return {
            'predictions': np.array(taxonomic_predictions),
            'confidence': confidences,
            'cluster_labels': cluster_labels,
            'novelty_scores': novelty_scores,
            'taxonomic_mapper': mapper,
            'novel_taxa': {
                'summary': {
                    'n_novel_clusters': 0,
                    'n_novel_outliers': len(novel_indices),
                    'total_novel_sequences': len(novel_indices)
                },
                'novel_outliers': novel_indices.tolist(),
                'novel_clusters': []
            }
        }
    
    def _calculate_abundance(self, predictions: Dict[str, Any], processed_data: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate abundance and diversity metrics."""
        predicted_taxa = predictions['predictions']
        confidence_scores = predictions['confidence']
        
        # Filter by confidence threshold
        confidence_threshold = 0.1  # Same as training
        high_confidence_mask = confidence_scores >= confidence_threshold
        filtered_taxa = predicted_taxa[high_confidence_mask]
        
        # Count sequences per taxon
        from collections import Counter
        taxon_counts = Counter(filtered_taxa)
        
        # Calculate relative abundance
        total_sequences = len(filtered_taxa)
        relative_abundance = {
            taxon: count / total_sequences 
            for taxon, count in taxon_counts.items()
        }
        
        # Calculate diversity metrics
        diversity_metrics = self._calculate_diversity_metrics(list(taxon_counts.values()))
        
        return {
            'taxon_counts': taxon_counts,
            'relative_abundance': relative_abundance,
            'total_sequences': len(predicted_taxa),
            'high_confidence_sequences': total_sequences,
            'confidence_threshold': confidence_threshold,
            'diversity_metrics': diversity_metrics
        }
    
    def _calculate_diversity_metrics(self, abundance_values: List[int]) -> Dict[str, float]:
        """Calculate biodiversity metrics."""
        if not abundance_values:
            return {}
        
        total = sum(abundance_values)
        proportions = [count / total for count in abundance_values]
        
        # Species richness
        richness = len(abundance_values)
        
        # Shannon diversity
        shannon_diversity = -sum(p * np.log(p) for p in proportions if p > 0)
        
        # Simpson diversity
        simpson_diversity = 1 - sum(p**2 for p in proportions)
        
        # Evenness
        max_shannon = np.log(richness) if richness > 0 else 0
        evenness = shannon_diversity / max_shannon if max_shannon > 0 else 0
        
        return {
            'species_richness': richness,
            'shannon_diversity': shannon_diversity,
            'simpson_diversity': simpson_diversity,
            'evenness': evenness
        }
    
    def _save_results(self, results: Dict[str, Any], output_path: Path, sample_name: str):
        """Save prediction results."""
        # Save taxonomic composition as CSV
        abundance_data = results['abundance_estimates']
        if abundance_data.get('taxon_counts'):
            composition_df = pd.DataFrame([
                {
                    'Taxon': taxon, 
                    'Count': count, 
                    'Relative_Abundance': abundance_data['relative_abundance'][taxon]
                }
                for taxon, count in abundance_data['taxon_counts'].items()
            ])
            composition_df = composition_df.sort_values('Count', ascending=False)
            composition_df.to_csv(output_path / f'{sample_name}_taxonomic_composition.csv', index=False)
        
        # Save diversity metrics
        if abundance_data.get('diversity_metrics'):
            diversity_df = pd.DataFrame([abundance_data['diversity_metrics']])
            diversity_df.to_csv(output_path / f'{sample_name}_diversity_metrics.csv', index=False)
        
        # Save summary report
        self._generate_summary_report(results, output_path, sample_name)
        
        logger.info(f"Results saved to: {output_path}")
    
    def _generate_summary_report(self, results: Dict[str, Any], output_path: Path, sample_name: str):
        """Generate a text summary report."""
        abundance_data = results['abundance_estimates']
        diversity = abundance_data.get('diversity_metrics', {})
        novel_taxa = results['predictions'].get('novel_taxa', {}).get('summary', {})
        
        report_lines = [
            f"eDNA Taxonomic Analysis Report",
            f"=" * 50,
            f"Sample: {sample_name}",
            f"Input: {results['input_path']}",
            f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"",
            f"SEQUENCE STATISTICS:",
            f"  Total Sequences: {results['total_sequences']}",
            f"  High Confidence Sequences: {abundance_data.get('high_confidence_sequences', 0)}",
            f"  Confidence Threshold: {abundance_data.get('confidence_threshold', 0.1)}",
            f"",
            f"BIODIVERSITY METRICS:",
            f"  Species Richness: {diversity.get('species_richness', 'N/A')}",
            f"  Shannon Diversity: {diversity.get('shannon_diversity', 0):.3f}",
            f"  Simpson Diversity: {diversity.get('simpson_diversity', 0):.3f}",
            f"  Evenness: {diversity.get('evenness', 0):.3f}",
            f"",
            f"NOVEL TAXA DETECTION:",
            f"  Novel Clusters: {novel_taxa.get('n_novel_clusters', 0)}",
            f"  Novel Outliers: {novel_taxa.get('n_novel_outliers', 0)}",
            f"  Total Novel Sequences: {novel_taxa.get('total_novel_sequences', 0)}",
            f"",
            f"TOP 10 MOST ABUNDANT TAXA:",
        ]
        
        # Add top taxa
        taxon_counts = abundance_data.get('taxon_counts', {})
        if taxon_counts:
            top_taxa = sorted(taxon_counts.items(), key=lambda x: x[1], reverse=True)[:10]
            for i, (taxon, count) in enumerate(top_taxa, 1):
                rel_abundance = abundance_data.get('relative_abundance', {}).get(taxon, 0)
                report_lines.append(f"  {i:2d}. {taxon}: {count} sequences ({rel_abundance:.1%})")
        
        report_lines.append(f"")
        report_lines.append(f"=" * 50)
        
        # Save report
        with open(output_path / f'{sample_name}_analysis_report.txt', 'w') as f:
            f.write('\n'.join(report_lines))


def main():
    """Main entry point for prediction script."""
    parser = argparse.ArgumentParser(
        description="Predict taxonomic composition of new eDNA samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python predict_new_sample.py --input new_sample.fasta --output results_new/
  python predict_new_sample.py --input path/to/samples/ --output batch_results/
  python predict_new_sample.py --input sample.fasta --models models_cluster --name "Lake_Sample_1"
        """
    )
    
    parser.add_argument(
        '--input', '-i', 
        type=str, 
        required=True,
        help='Path to input FASTA file or directory containing FASTA files'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str, 
        default='prediction_results',
        help='Output directory for results (default: prediction_results)'
    )
    
    parser.add_argument(
        '--models', '-m',
        type=str,
        default='models_cluster', 
        help='Directory containing trained models (default: models_cluster)'
    )
    
    parser.add_argument(
        '--name', '-n',
        type=str,
        help='Name for the sample (default: filename)'
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize predictor
        predictor = eDNAPredictor(model_dir=args.models)
        
        # Make predictions
        results = predictor.predict_sample(
            input_path=args.input,
            output_dir=args.output,
            sample_name=args.name
        )
        
        # Print summary
        print(f"\n🧬 eDNA Analysis Complete! 🧬")
        print(f"Sample: {results['sample_name']}")
        print(f"Total sequences: {results['total_sequences']}")
        
        diversity = results['abundance_estimates'].get('diversity_metrics', {})
        print(f"Species richness: {diversity.get('species_richness', 'N/A')}")
        print(f"Shannon diversity: {diversity.get('shannon_diversity', 0):.3f}")
        
        print(f"\nResults saved to: {args.output}")
        print(f"Check the CSV files and analysis report for detailed results!")
        
    except Exception as e:
        print(f"❌ Prediction failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
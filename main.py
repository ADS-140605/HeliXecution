#!/usr/bin/env python3
"""
Main execution script for the eDNA AI Pipeline.
Integrates all components for taxonomic classification and biodiversity assessment.
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path
from typing import Dict, Any
import warnings
warnings.filterwarnings('ignore')

# Add src to path
sys.path.append(str(Path(__file__).parent / 'src'))

# Import pipeline components
from src.models.ensemble_pipeline import eDNAEnsemblePipeline
from src.visualization.biodiversity_plots import BiodiversityVisualizer

def setup_logging(config: Dict[str, Any]):
    """
    Set up logging configuration.
    
    Args:
        config: Configuration dictionary
    """
    log_config = config.get('logging', {})
    log_level = getattr(logging, log_config.get('level', 'INFO'))
    log_file = log_config.get('log_file', 'logs/pipeline.log')
    console_output = log_config.get('console_output', True)
    
    # Create logs directory if it doesn't exist
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Configure logging
    handlers = [logging.FileHandler(log_file)]
    if console_output:
        handlers.append(logging.StreamHandler(sys.stdout))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    
    # Set specific loggers to WARNING to reduce noise
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('plotly').setLevel(logging.WARNING)

def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Configuration dictionary
    """
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        print(f"Error loading configuration: {e}")
        sys.exit(1)

def create_output_directories(config: Dict[str, Any]):
    """
    Create output directories based on configuration.
    
    Args:
        config: Configuration dictionary
    """
    output_config = config.get('output', {})
    
    directories = [
        output_config.get('model_save_dir', 'models'),
        output_config.get('results_dir', 'results'),
        output_config.get('plots_dir', 'plots'),
        output_config.get('reports_dir', 'reports'),
        'logs'
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)

def run_pipeline(data_path: str, config: Dict[str, Any], blast_db_path: str = None):
    """
    Run the complete eDNA analysis pipeline.
    
    Args:
        data_path: Path to input data (FASTA files)
        config: Configuration dictionary
        blast_db_path: Optional path to BLAST database
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting eDNA AI Pipeline")
    logger.info(f"Input data: {data_path}")
    logger.info(f"BLAST database: {blast_db_path if blast_db_path else 'None'}")
    
    # Initialize pipeline
    model_config = config.get('models', {})
    ensemble_config = config.get('ensemble', {})
    
    pipeline = eDNAEnsemblePipeline(
        use_cnn=model_config.get('cnn', {}).get('enabled', True),
        use_gnn=model_config.get('gnn', {}).get('enabled', True),
        use_rf=model_config.get('random_forest', {}).get('enabled', True),
        use_clustering=model_config.get('clustering', {}).get('enabled', True),
        ensemble_method=ensemble_config.get('method', 'weighted_voting'),
        confidence_threshold=ensemble_config.get('confidence_threshold', 0.7)
    )
    
    # Set model weights
    model_weights = ensemble_config.get('model_weights', {})
    if model_weights:
        pipeline.model_weights.update(model_weights)
    
    try:
        # Step 1: Load and preprocess data
        logger.info("Step 1: Loading and preprocessing data")
        data_proc_config = config.get('data_processing', {})
        
        processed_data = pipeline.load_and_preprocess_data(
            data_path=data_path,
            blast_db_path=blast_db_path,
            sample_limit=data_proc_config.get('sample_limit')
        )
        
        logger.info(f"Loaded {len(processed_data['sequences'])} sequences")
        
        # Step 2: Train models
        logger.info("Step 2: Training models")
        training_config = config.get('training', {})
        
        training_results = pipeline.train_models(
            processed_data=processed_data,
            test_size=training_config.get('test_size', 0.2),
            validation_size=training_config.get('validation_size', 0.1),
            **model_config  # Pass model-specific configurations
        )
        
        # Step 3: Make predictions
        logger.info("Step 3: Making predictions")
        predictions = pipeline.predict(processed_data)
        
        # Step 4: Estimate abundance and diversity
        logger.info("Step 4: Estimating abundance and diversity")
        abundance_estimates = pipeline.estimate_abundance(processed_data, predictions)
        
        # Step 5: Generate visualizations and reports
        logger.info("Step 5: Generating visualizations and reports")
        
        # Initialize visualizer
        viz_config = config.get('visualization', {})
        visualizer = BiodiversityVisualizer(
            style=viz_config.get('style', 'seaborn'),
            figsize=tuple(viz_config.get('figsize', [12, 8]))
        )
        
        # Create comprehensive results dictionary
        pipeline_results = {
            'pipeline_configuration': {
                'models_used': {
                    'cnn': pipeline.use_cnn,
                    'gnn': pipeline.use_gnn,
                    'random_forest': pipeline.use_rf,
                    'clustering': pipeline.use_clustering
                },
                'ensemble_method': pipeline.ensemble_method,
                'confidence_threshold': pipeline.confidence_threshold
            },
            'training_results': training_results,
            'predictions': predictions,
            'abundance_estimates': abundance_estimates
        }
        
        # Generate comprehensive report
        output_config = config.get('output', {})
        reports_dir = output_config.get('reports_dir', 'reports')
        
        if viz_config.get('generate_html_report', True):
            report_path = visualizer.generate_comprehensive_report(
                pipeline_results=pipeline_results,
                processed_data=processed_data,
                output_dir=reports_dir
            )
            logger.info(f"Comprehensive report generated: {report_path}")
        
        # Save models
        if output_config.get('save_models', True):
            model_save_dir = output_config.get('model_save_dir', 'models')
            pipeline.save_models(model_save_dir)
            logger.info(f"Models saved to: {model_save_dir}")
        
        # Save results
        results_dir = output_config.get('results_dir', 'results')
        results_path = Path(results_dir)
        results_path.mkdir(parents=True, exist_ok=True)
        
        # Save key results as CSV files
        if abundance_estimates.get('taxon_counts'):
            import pandas as pd
            
            # Taxonomic composition
            composition_df = pd.DataFrame([
                {'Taxon': taxon, 'Count': count, 'Relative_Abundance': abundance_estimates['relative_abundance'][taxon]}
                for taxon, count in abundance_estimates['taxon_counts'].items()
            ])
            composition_df.to_csv(results_path / 'taxonomic_composition.csv', index=False)
            
            # Diversity metrics
            diversity_df = pd.DataFrame([abundance_estimates['diversity_metrics']])
            diversity_df.to_csv(results_path / 'diversity_metrics.csv', index=False)
            
            logger.info(f"Results saved to: {results_dir}")
        
        # Print summary statistics
        print_summary(abundance_estimates, predictions)
        
        logger.info("eDNA AI Pipeline completed successfully!")
        
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}")
        raise

def print_summary(abundance_estimates: Dict[str, Any], predictions: Dict[str, Any]):
    """
    Print summary statistics to console.
    
    Args:
        abundance_estimates: Abundance estimation results
        predictions: Prediction results
    """
    print("\n" + "="*60)
    print("                eDNA ANALYSIS SUMMARY")
    print("="*60)
    
    # Basic statistics
    total_sequences = abundance_estimates.get('total_sequences', 0)
    high_conf_sequences = abundance_estimates.get('high_confidence_sequences', 0)
    
    print(f"Total Sequences Processed: {total_sequences}")
    print(f"High Confidence Sequences: {high_conf_sequences}")
    print(f"Confidence Threshold: {abundance_estimates.get('confidence_threshold', 'N/A')}")
    
    # Diversity metrics
    diversity = abundance_estimates.get('diversity_metrics', {})
    print(f"\nBIODIVERSITY METRICS:")
    print(f"  Species Richness: {diversity.get('species_richness', 'N/A')}")
    shannon = f"{diversity['shannon_diversity']:.3f}" if 'shannon_diversity' in diversity else 'N/A'
    simpson = f"{diversity['simpson_diversity']:.3f}" if 'simpson_diversity' in diversity else 'N/A'
    evenness = f"{diversity['evenness']:.3f}" if 'evenness' in diversity else 'N/A'
    print(f"  Shannon Diversity: {shannon}")
    print(f"  Simpson Diversity: {simpson}")
    print(f"  Evenness: {evenness}")
    
    # Novel taxa
    novel_taxa = predictions.get('novel_taxa', {}).get('summary', {})
    print(f"\nNOVEL TAXA DETECTION:")
    print(f"  Novel Clusters: {novel_taxa.get('n_novel_clusters', 0)}")
    print(f"  Novel Outliers: {novel_taxa.get('n_novel_outliers', 0)}")
    print(f"  Total Novel Sequences: {novel_taxa.get('total_novel_sequences', 0)}")
    
    # Top taxa
    taxon_counts = abundance_estimates.get('taxon_counts', {})
    if taxon_counts:
        print(f"\nTOP 5 MOST ABUNDANT TAXA:")
        top_taxa = sorted(taxon_counts.items(), key=lambda x: x[1], reverse=True)[:5]
        for i, (taxon, count) in enumerate(top_taxa, 1):
            rel_abundance = abundance_estimates.get('relative_abundance', {}).get(taxon, 0)
            print(f"  {i}. {taxon}: {count} sequences ({rel_abundance:.1%})")
    
    print("="*60)

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="eDNA AI Pipeline for Taxonomic Classification and Biodiversity Assessment",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --data /path/to/sequences.fasta
  python main.py --data /path/to/fasta_directory --config custom_config.yaml
  python main.py --data /path/to/sequences.fasta --blast-db /path/to/blastdb --config config.yaml
        """
    )
    
    parser.add_argument(
        '--data', 
        type=str, 
        required=True,
        help='Path to input FASTA file or directory containing FASTA files'
    )
    
    parser.add_argument(
        '--config', 
        type=str, 
        default='config/pipeline_config.yaml',
        help='Path to configuration YAML file (default: config/pipeline_config.yaml)'
    )
    
    parser.add_argument(
        '--blast-db', 
        type=str, 
        default=None,
        help='Path to BLAST database directory (optional)'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='eDNA AI Pipeline v1.0.0'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Set up logging
    setup_logging(config)
    
    # Create output directories
    create_output_directories(config)
    
    # Run pipeline
    try:
        run_pipeline(
            data_path=args.data,
            config=config,
            blast_db_path=args.blast_db
        )
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nPipeline failed: {str(e)}")
        logging.error(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
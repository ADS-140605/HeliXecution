# 🧬 How to Predict New eDNA Samples

Your eDNA AI pipeline is now trained and ready to analyze new samples! Here's how to use it:

## 📁 Quick Start

### 1. **Single Sample Analysis**
```bash
python predict_new_sample.py --input your_sample.fasta --output results/ --name "Lake_Sample_1"
```

### 2. **Batch Analysis (Multiple Samples)**
```bash
python predict_new_sample.py --input samples_folder/ --output batch_results/
```

### 3. **Using Different Model Directory**
```bash
python predict_new_sample.py --input sample.fasta --models models_blast --output results/
```

## 📋 Input Requirements

### **FASTA File Format**
Your new eDNA sample should be in FASTA format:
```
>Sequence_1
ATCGATCGATCGATCGATCGATCGATCG...
>Sequence_2  
GCTAGCTAGCTAGCTAGCTAGCTAGCT...
>Sequence_3
AAATTTGGCCAAATTTGGCCAAATTT...
```

### **File Requirements**
- **Format**: FASTA (.fasta, .fa, .fas)
- **Sequence Length**: 100-2000 nucleotides (same as training)
- **Quality**: <5% ambiguous nucleotides (N's)
- **Size**: Any number of sequences (from 1 to thousands)

## 📊 Output Files

For each sample, you'll get:

### **1. Taxonomic Composition** (`SampleName_taxonomic_composition.csv`)
```csv
Taxon,Count,Relative_Abundance
Mixed_Escherichia_Gemmata,45,0.225
Mixed_Anabaena_Pseudomonas,38,0.19
Bartonella,25,0.125
...
```

### **2. Diversity Metrics** (`SampleName_diversity_metrics.csv`)
```csv
species_richness,shannon_diversity,simpson_diversity,evenness
12,2.485,0.923,0.998
```

### **3. Analysis Report** (`SampleName_analysis_report.txt`)
```
eDNA Taxonomic Analysis Report
==================================================
Sample: Lake_Sample_1
Total Sequences: 200
Species Richness: 12
Shannon Diversity: 2.485

TOP 10 MOST ABUNDANT TAXA:
  1. Mixed_Escherichia_Gemmata: 45 sequences (22.5%)
  2. Mixed_Anabaena_Pseudomonas: 38 sequences (19.0%)
  3. Bartonella: 25 sequences (12.5%)
  ...
```

## 🎯 Command Options

| Option | Short | Description | Example |
|--------|-------|-------------|---------|
| `--input` | `-i` | Input FASTA file or directory | `--input sample.fasta` |
| `--output` | `-o` | Output directory | `--output results/` |
| `--name` | `-n` | Custom sample name | `--name "River_Site_A"` |
| `--models` | `-m` | Model directory | `--models models_cluster` |

## 📖 Real Examples

### **Example 1: Lake Water Sample**
```bash
# Analyze lake water eDNA
python predict_new_sample.py \
  --input lake_water_june2024.fasta \
  --output lake_analysis/ \
  --name "Lake_June_2024"
```

**Expected Results:**
- Identifies bacterial communities like Cyanobacteria, Proteobacteria
- Calculates diversity metrics for the lake ecosystem  
- Detects potential novel taxa unique to the lake

### **Example 2: Soil Sample**
```bash
# Analyze soil eDNA
python predict_new_sample.py \
  --input forest_soil.fasta \
  --output soil_results/ \
  --name "Forest_Soil_Sample"
```

**Expected Results:**
- Different bacterial profile (more Actinobacteria, Acidobacteria)
- Higher diversity typical of soil environments
- Detection of soil-specific bacterial genera

### **Example 3: Marine Sample**
```bash
# Analyze marine eDNA
python predict_new_sample.py \
  --input ocean_sample.fasta \
  --output marine_analysis/ \
  --name "Deep_Sea_Sample"
```

**Expected Results:**
- Marine-specific bacteria (Vibrio, Pseudomonas, etc.)
- Novel taxa detection for unexplored deep-sea environments
- Marine biodiversity metrics

## 🚀 Advanced Usage

### **Batch Processing Multiple Samples**
```bash
# Directory structure:
samples/
├── lake_1.fasta
├── lake_2.fasta
├── soil_1.fasta
└── marine_1.fasta

# Process all at once:
python predict_new_sample.py --input samples/ --output batch_results/
```

### **Comparison Analysis**
```bash
# Compare different environments
python predict_new_sample.py --input lake.fasta --output comparison/ --name "Lake"
python predict_new_sample.py --input soil.fasta --output comparison/ --name "Soil"
python predict_new_sample.py --input marine.fasta --output comparison/ --name "Marine"

# Then compare the CSV files to see differences in:
# - Species richness
# - Shannon diversity
# - Dominant taxa
# - Novel taxa detection
```

## 📈 Interpreting Results

### **Taxonomic Names**
- **Single names** (e.g., "Bartonella"): High confidence assignment
- **Mixed names** (e.g., "Mixed_Escherichia_Gemmata"): Cluster contains multiple similar genera
- **Novel detection**: Sequences that don't match known patterns

### **Diversity Metrics**
- **Species Richness**: Number of different taxa (higher = more diverse)
- **Shannon Diversity**: Accounts for evenness (0-5, higher = more diverse)
- **Simpson Diversity**: Probability two sequences are different (0-1, higher = more diverse)
- **Evenness**: How evenly distributed taxa are (0-1, higher = more even)

### **Confidence Interpretation**
- **High confidence sequences**: Match well to known patterns
- **Novel sequences**: Potentially new or rare species
- **Low confidence**: Uncertain assignments (still reported)

## 🔧 Troubleshooting

### **Common Issues**

**❌ "No FASTA files found"**
- Check file extensions (.fasta, .fa, .fas)
- Verify file path is correct

**❌ "All sequences filtered out"**
- Sequences too short (<100 bp) or too long (>2000 bp)
- Too many ambiguous nucleotides (N's)

**❌ "Model not found"**
- Ensure `models_cluster/` directory exists
- Run the training pipeline first

### **Sample Size Guidelines**
- **Small samples** (1-50 sequences): Basic analysis, fewer clusters
- **Medium samples** (50-500 sequences): Good diversity analysis
- **Large samples** (500+ sequences): Full taxonomic profiling, reliable statistics

## 🎯 What Makes This Different

Unlike traditional BLAST-based approaches, your AI pipeline:

✅ **Identifies novel taxa** not in databases
✅ **Groups similar sequences** into meaningful clusters  
✅ **Provides bacterial genus names** instead of just "OTU_1, OTU_2"
✅ **Calculates diversity metrics** automatically
✅ **Works with any eDNA sample** from any environment

## 📞 Next Steps

1. **Try it with your real samples!** 
2. **Compare different environments** (lake vs. soil vs. marine)
3. **Track changes over time** (seasonal variation)
4. **Identify novel taxa** for further investigation

Your trained model is ready to discover the microbial world in your eDNA samples! 🦠✨
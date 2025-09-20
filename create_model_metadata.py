#!/usr/bin/env python3
"""
Create model metadata for prediction system
"""

import joblib
from pathlib import Path

# Create metadata for our clustering model
metadata = {
    'use_cnn': False,
    'use_gnn': False, 
    'use_rf': False,
    'use_clustering': True,
    'ensemble_method': 'weighted_voting',
    'confidence_threshold': 0.1,
    'model_weights': {
        'cnn': 0.0,
        'gnn': 0.0,
        'rf': 0.0
    }
}

# Save to models_cluster directory
models_dir = Path("models_cluster")
models_dir.mkdir(exist_ok=True)

joblib.dump(metadata, models_dir / "pipeline_metadata.joblib")
print("Model metadata saved to models_cluster/pipeline_metadata.joblib")
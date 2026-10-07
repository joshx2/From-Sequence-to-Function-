# Protein-prediction
interactive AI-based system for predicting protein function from protein sequence.
Senior design project: predict protein function (GO terms) from sequences using ESM-2 embeddings, served through an interactive dashboard.

## Repo layout

- 'data/' - raw and processed datasets (gitignored, never commit)
- 'src/' - preprocessing, model, training, evaluation code
- 'api/' - FastAPI backend that serves predictions
- 'dashboard/' - frontend (HTML/JS/Plotly)
- 'notebooks/' - exploration and one-off analysis
- 'scripts/' - setup/verification utilities

## Setup

'''bash
conda create -n protein-func python=3.10 -y
conda activate protein-func
pip install -r requirements.txt
'''

then verify everything works:

'''bash
python scripts/smoke_test.py
'''

## Data sources

- CAFA6 (Kaggle, recommended starting point - already bundled sequences + GO terms):  https://www.kaggle.com/competitions/cafa-6-protein-function-prediction/data?utm_source=chatgpt.com
- Swiss-Prot flat file (more parsing work, but includes everything current):  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
- GO ontology (needed for label propagation with goatools):
  https://purl.obolibrary.org/obo/go/go-basic.obo

## Status

Repo scaffold only - pipeline code comes next as we work through the data prep, embedding extraction, model training, and dashboard steps. 

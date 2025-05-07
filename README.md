# Probabilistic-Genome-Alignment

This project implements a custom algorithm for aligning a short query sequence to a long **probabilistic genome**, where each nucleotide in the database has an associated confidence score. Unlike BLAST, which is designed for deterministic genomes, this algorithm incorporates the uncertainty in nucleotide predictions and uses a modified Needleman-Wunsch extension to produce a single high-likelihood alignment.

## Project Overview

The alignment process is broken into five main steps:

1. **Data Preprocessing**: Loads the database and its associated probabilities, then generates word-index dictionaries.
2. **Subsequence Identification**: Identifies promising subarrays using hit counts and uniqueness rewards.
3. **Targeted Reprocessing**: Rebuilds word dictionaries for the selected subarray using probabilistic recombination.
4. **Diagonal Filtering**: Narrows down the alignment region using hit density across diagonals.
5. **Probabilistic Alignment**: Aligns the query using a log-likelihood scoring version of Needleman-Wunsch.

The algorithm is optimized for runtime and accuracy using heuristics, threading and empirical parameter tuning.

### Folder Structure

```
Probabilistic-Genome-Alignment/
├── README.md                 ← This file
├── LICENSE                   ← All rights reserved license
├── report/
│   └── Mélanie Ghaby Final Report.pdf   ← Full project report
├── data/
│   ├── predicted_database.txt
│   └── predicted_database_probs.txt
├── results/
│   ├── analyse.R
│   └── official_dataset.xlsx
├── src/
│   ├── alignment.py
│   ├── subarray_identifier.py
│   ├── diagonal_dictionary_recomposer.py
│   ├── targeted_reprocessing.py
│   ├── data_preprocessing.py
│   ├── code.py
│   └── query_generator_and_tester.py
├── json_data/
│   ├── w_6.json
│   ├── w_7.json
│   └── ...
```

## Requirements

- Python 3.9+
- Packages:
  - `numpy`
  - `pandas`
  - `openpyxl`
  - `concurrent.futures` (standard in Python 3.9)
- (Optional) R 4.4.1+ for analyzing simulation results

## How to Run

1. Place your query and database files in the `data/` folder.
2. Edit and run `src/code.py` to execute the full pipeline.
3. Output Excel results will appear in the `results/` folder.

## Report

For full methodology, heuristics, and evaluation results, see:
`report/Mélanie Ghaby Final Report.pdf`

## License

**Copyright © 2025 Mélanie Ghaby. All rights reserved.**

This code and accompanying materials may **not** be copied, distributed, modified, or used in any form without prior written permission from the author.  
This project is shared publicly for **viewing and demonstration purposes only**.

Please replace the contents of this file with relevant instructions for your repository or remove this file entirely.

This directory would generally contain source code files that contain the core code to implement the method and various utility/auxiliary functions.

Scripts/code files that execute the overall workflow to carry out an analysis and generate results for the manuscript might be placed in the main directory.






## STRUCTURE runs (brook trout)

### Software
- STRUCTURE (version: <fill in>)
- Binary location: ./structure  (or on PATH)

### Inputs
- Parameter files: code/structure/mainparams, code/structure/extraparams
- Genotype file: data/raw/brooktrout/brooktrout.txt
  (Note: mainparams must point to this file, or add -i if your STRUCTURE build requires it.)

### Outputs
- Written to: data/processed/structure_run/
- Naming: results_K{K}_rep{rep} (plus any STRUCTURE suffixes)

### Command used
```bash
mkdir -p data/processed/structure_run

for K in {1..10}; do
  for rep in {1..20}; do
    SEED=$((1000 * K + rep))
    ./structure \
      -K $K \
      -D $SEED \
      -m code/structure/mainparams \
      -e code/structure/extraparams \
      -o data/processed/structure_run/results_K${K}_rep${rep}
  done
done

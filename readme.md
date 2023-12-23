# Code & Data Availability for ___Group size modulates kinship dynamics and selection on social traits___

---

### contents

- [model](#model)
- [data](#data)
- [plot](#plot)
- [e.sh](#e.sh)

---

### model

This folder contains the Python script for calculating individuals' local genetic relatedness to others and evaluating the selective pressures on their social behaviours expressed in such locally dynamic social environments. The script `KINSHIP_DYNAMICS_AND_SELECTION.py` is used to generate the predictions presented in the manuscript (i.e., both kinship dynamics and selective pressures), while the script `KINSHIP_DYNAMICS.py` is for predicting the patterns of kinship dynamics only (if so desired). All results from the model are written to the _res_ folder by default.

---

### data

The killer whale data used in the manuscript (which will be made available shortly).

---

### plot

This folder contains the R script visualizing the results from both model predictions and empirical observations (main figures for the manuscript are hosted in `main` while supplementary figures are in `supp`).

---

### e.sh

The shell script to execute the scripts.
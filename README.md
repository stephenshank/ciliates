#Ciliates

This repository contains codes for analyzing telemore end-binding proteins found in ciliates.

After cloning this repository, make sure to obtain a copy of `ciliates_dataset.tar.gz` and place it in the `data` directory. Then run

```sh
bash pipeline.sh
```

Note that you may need to install requirements (for instance, using `pip install -r requirements.txt`). To clean, run

```sh
bash clean.sh
```

## Directory Structure

`data` - Contains all data for the analysis. Subdirectories are described below.

`data/initial` - Contains initial dataset, consisting of protein/nucleotide sequences and subunit memberships.
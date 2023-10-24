# Arizona searches by imputation

Code accompanying the manuscript *"Performing Arizona searches by imputation"* by Lappo & Rosenberg (2023)

## How to run

**0.** Make sure your machine has [`git`](https://git-scm.com/downloads) and [`git-lfs`](https://git-lfs.com) installed.

Then, clone this repository to your machine by running either of these commands.
```
gh repo clone EgorLappo/arizona-searches-by-imputation
git clone https://github.com/EgorLappo/arizona-searches-by-imputation.git
```

Make sure the input data has been downloaded by running `git lfs fetch`.

**1.** Install `nix` with [this installer](https://zero-to-nix.com/start/install). TLDR: run the following command

```bash
curl --proto '=https' --tlsv1.2 -sSf -L https://install.determinate.systems/nix | sh -s -- install
```

**2.** Run `nix develop` to enter a shell with all required dependencies and scripts. Run the `run-analysis` command in this shell to generate figures/tables. See `flake.nix` for more details. The `run-analysis` command will execute the following commands in order:

- `prepare` will process and index the input `.vcf` files, download the BEAGLE `.jar` file into the `beagle` folder.
- `impute` will generate 100 replicate train-test splits (`processed_data` folder), run BEAGLE's imputation algorithm, and summarize the results with `bcftools` (`imputation_results` folder).
- `matches` will read in the genotype data, compute all the required numbers of matches as well as the necessary theoretical quantities (`computed_matches` folder).
- finally, the R script `run_tests.R` will be run to perform the Wilcoxon tests and the python script `make_figures.py` will plot the figures (`make_figures` and `figures` folders).

Each of these steps can be performed manually by running the individual commands or scripts in the shell provided by `nix`.

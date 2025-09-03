# Pansim
A fast multi-threaded pangenome dynamics simulator written in Rust.

Pansim is a Wright-Fisher simulator, incorporating gene gain/loss, recombination and selection.

Pansim generates pairwise core vs. accessory distances, like [PopPUNK](https://github.com/bacpop/PopPUNK), core genome alignments and gene presence/absence matrices.

## Installation

First install rust (>=1.89.0). This can be done using [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html):

```
micromamba install conda-forge::rust
```

Then clone this repository and install pansim:

```
git clone https://github.com/bacpop/Pansim.git && cd Pansim/pansim && cargo build --release && cd ../..
```

The executable will be in `Pansim/pansim/target/release/pansim`. You can add this to your path to allow running at anytime:

```
export PATH="Pansim/pansim/target/release:$PATH"
```

To make this permanent, add the above line to your `~/.bashrc` or `~/.zshrc` file.

## Running Pansim

Pansim can be run by calling the executable in `Pansim/pansim/target/release/pansim`. If added to your path, simply run:

```
pansim --outpref test_run
```

All command line arguments can be found below:
```
USAGE:
    pansim [OPTIONS]

OPTIONS:
        --HGT_rate <HGT_rate>
            HGT rate, as number of accessory sites transferred per core genome mutation. [default:
            0.05]

        --HR_rate <HR_rate>
            Homologous recombination rate, as number of core sites transferred per core genome
            mutation. [default: 0.05]

        --avg_gene_freq <avg_gene_freq>
            Average proportion of genes in pangenome present in an individual. Includes core and
            accessory genes. [default: 0.5]

        --competition_strength <competition_strength>
            Strength of competition felt by strain to all others. 0.0 = no competition [default:
            0.0]

        --core_genes <core_genes>
            Number of core genes in pangenome. [default: 2000]

        --core_mu <core_mu>
            Average core SNP mutation rate (per site per genome per generation in core genome). Must
            be > 0.0. [default: 0.05]

        --core_size <core_size>
            Number of nucleotides in core genome. [default: 1200000]

        --genome_size_penalty <genome_size_penalty>
            Multiplier for each gene difference between avg_gene_freq and observed value. [default:
            0.99]

    -h, --help
            Print help information

        --max_distances <max_distances>
            Maximum number of pairwise distances to calculate. [default: 100000]

        --n_gen <n_gen>
            Number of generations to simulate. [default: 100]

        --neg_lambda <neg_lambda>
            Lambda value for exponential distribution of negatively selected genes. Must be > 0.0.
            [default: 10.0]

        --no_control_genome_size
            Removes penalisation of genome sizes deviating from average.

        --outpref <outpref>
            Output prefix path. [default: distances]

        --pan_genes <pan_genes>
            Total number of genes in pangenome (core + accessory). [default: 6000]

        --pop_size <pop_size>
            Number of individuals in population. [default: 1000]

        --pos_lambda <pos_lambda>
            Lambda value for exponential distribution of positively selected genes. Must be > 0.0.
            [default: 10.0]

        --print_dist
            Print per-generation average pairwise distances.

        --print_matrices
            Prints core and accessory matrices.

        --print_selection
            Prints selection coefficients.

        --prop_genes2 <prop_genes2>
            Proportion of pangenome made up of compartment 2 genes. Must be 0.0 <= X <= 1.0.
            [default: 0.1]

        --prop_positive <prop_positive>
            Proportion of pangenome made up of positively selected genes. Must be 0.0 <= X <= 1.0.
            If negative, neutral selection is simulated. [default: -0.1]

        --rate_genes1 <rate_genes1>
            Average number of accessory genes that are gained/lost per site per genome per
            generation in gene compartment 1. Must be >= 0.0. [default: 1.0]

        --rate_genes2 <rate_genes2>
            Average number of accessory genes that are gained/lost per site per genome per
            generation in gene compartment 2. Must be >= 0.0. [default: 1000.0]

        --seed <seed>
            Seed for random number generation. [default: 0]

        --threads <threads>
            Number of threads. [default: 1]

    -V, --version
            Print version information

        --verbose
            Prints per generation information.
```

## Output


### Core vs. accessory distances

By default, Pansim generates core and accessory pairwise distance, named `<outpref>.tsv`, like [PopPUNK](https://github.com/bacpop/PopPUNK). The core and accessory distances per comparison are the first and second column respectively.

### Gene frequencies

By default, Pansim generates a gene frequencies file, named `<outpref>_freqs.txt`, where each line designates a gene in the pangenome and the proportions of total individuals it is found in.

### Selection coefficients

Specifying `--print_selection` generates `<outpref>_selection.tsv`, which is a distribution of selection coefficients for all genes in the simulation.

### Core genome alignment and gene presenence/absence matrix

Specifying `--print_matrices` generates `<outpref>_core_genome.csv` and `<outpref>_pangenome.csv`, which are the core genome alignment and gene presence/absence matrix respectively. Each row is a single individual, each column is a site in the core genome/ gene in the pangenome.

### Per generation genome distances

Specifying `--print_dist` generates `<outpref>_per_gen.tsv`, which decribes the changes in core and accessory genome pairwise distances across the simulation time. Each row is a single generation. The columns denote (moving left to right): average core distance, core distance standard deviation, average accessory distance and accessory distance standard deviation.
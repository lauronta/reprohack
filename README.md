# Reprohackathon 2024
L'objectif de ce projet est de reproduire les analyses réalisées dans la publication suivante en utilisant les technologies de conteneurisation [Singularity](https://docs.sylabs.io/guides/latest/user-guide/#) et de gestion de pipeline [Snakemake](https://snakemake.github.io/). 
>  Peyrusson, F., Varet, H., Nguyen, T.K. et al. Intracellular Staphylococcus aureus persisters upon antibiotic exposure. Nat Commun 11, 2200 (2020). [https://doi.org/10.1038/s41467-020-15966-7](https://doi.org/10.1038/s41467-020-15966-7)

## Installation


### Prérequis :

- `singularity >= 3.0.1`
- `snakemake`

### Cloner le dépôt

    git clone https://github.com/lauronta/reprohack.git

### Lancer le workflow

    cd reprohack/
    ./run.sh
    
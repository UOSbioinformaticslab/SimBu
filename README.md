# SimBu Project

## Overview
SimBu is an R package designed for simulating single-cell RNA sequencing data and evaluating various batch correction methods. The package provides tools for generating synthetic datasets, introducing batch effects, and assessing the performance of different correction techniques.

## Installation
To install the SimBu package, you can use the following command in R:

```R
# Install the devtools package if you haven't already
install.packages("devtools")

# Install SimBu from GitHub
devtools::install_github("stef1949/SimBu")
```

## Usage
After installing the package, you can load it and use the functions provided for data generation and batch correction evaluation. Here is a simple example:

```R
library(SimBu)

# Generate a base single-cell dataset
base_dataset <- create_base_sc_dataset(n_genes = 1000, n_cells = 300)

# Generate simulations
simulation_list <- generate_simulations(simbu_dataset = base_dataset, n_sims = 3, n_samples_per_sim = 30)

# Introduce batch effects
batch_effect_data <- introduce_batch_effect(simulation_list)

# Evaluate batch correction methods
# (Add your evaluation code here)
```

## Contribution
Contributions are welcome! If you would like to contribute to the SimBu project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Make your changes and commit them.
4. Push your branch to your forked repository.
5. Create a pull request describing your changes.

Please ensure that your code adheres to the project's coding standards and includes appropriate tests.

## License
This project is licensed under the MIT License. See the LICENSE file for more details.
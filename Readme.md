# Python Code for calculating liklihood of observed MSA given a phylogenetic tree

This repository contains Python code for calculating likelihood (lkl) for phylogenetic trees. The code includes two main functionalities:

1. **Get LKL for All Nodes:**
    ```bash
    python3 run_lkl_calculation.py \
        --tree_path data/table.dat \
        --branch_lengths_path data/branchlength.dat \
        --msa_path data/msa.dat \
        --mu 1
    ```

2. **Get LKL for Root Node Only:**
    ```bash
    python3 run_lkl_calculation.py \
        --tree_path data/table.dat \
        --branch_lengths_path data/branchlength.dat \
        --msa_path data/msa.dat \
        --mu 1 \
        --output_lkl_root_node_only
    ```

### Test Datasets

To run the test datasets, you can use the following commands:

```bash
# Install required dependencies
pip install tabulate

# Run the test script
python3 test.py
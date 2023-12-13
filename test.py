import subprocess
from tabulate import tabulate

datasets = [
    ("msa.dat", "branchlength.dat", "table.dat", "1.09", "0.3", "-180.4314"),
    ("msa.dat", "branchlength.dat", "table.dat", "1.09", "1", "-182.0558"),
    ("msa.dat", "branchlength.dat", "table.dat", "1.09", "2", "-197.3004"),
    ("msa.dat", "branchlength.dat", "table.dat", "1.09", "0.5321", "-176.5073"),
    ("ENSG00000112282_MED23_NT.msa.dat", "ENSG00000112282_MED23_NT.branchlength.dat", "ENSG00000112282_MED23_NT.table.dat", "3.044", "0.3", "-11350.558"),
    ("ENSG00000112282_MED23_NT.msa.dat", "ENSG00000112282_MED23_NT.branchlength.dat", "ENSG00000112282_MED23_NT.table.dat", "3.044", "1", "-10229.941"),
    ("ENSG00000112282_MED23_NT.msa.dat", "ENSG00000112282_MED23_NT.branchlength.dat", "ENSG00000112282_MED23_NT.table.dat", "3.044", "2", "-10534.2058"),
    ("ENSG00000112282_MED23_NT.msa.dat", "ENSG00000112282_MED23_NT.branchlength.dat", "ENSG00000112282_MED23_NT.table.dat", "3.044", "1.1549", "-10211.8688"),
    ("ENSG00000112984_KIF20A_NT.msa.dat", "ENSG00000112984_KIF20A_NT.branchlength.dat", "ENSG00000112984_KIF20A_NT.table.dat", "3.525", "0.3", "-7412.2796"),
    ("ENSG00000112984_KIF20A_NT.msa.dat", "ENSG00000112984_KIF20A_NT.branchlength.dat", "ENSG00000112984_KIF20A_NT.table.dat", "3.525", "1", "-7134.6171"),
    ("ENSG00000112984_KIF20A_NT.msa.dat", "ENSG00000112984_KIF20A_NT.branchlength.dat", "ENSG00000112984_KIF20A_NT.table.dat", "3.525", "2", "-8029.2966"),
    ("ENSG00000112984_KIF20A_NT.msa.dat", "ENSG00000112984_KIF20A_NT.branchlength.dat", "ENSG00000112984_KIF20A_NT.table.dat", "3.525", "0.7021", "-7043.9356"),
    ("ENSG00000013016_EHD3_NT.msa.dat", "ENSG00000013016_EHD3_NT.branchlength.dat", "ENSG00000013016_EHD3_NT.table.dat", "5.209", "0.3", "-6802.5391"),
    ("ENSG00000013016_EHD3_NT.msa.dat", "ENSG00000013016_EHD3_NT.branchlength.dat", "ENSG00000013016_EHD3_NT.table.dat", "5.209", "1", "-7324.2722"),
    ("ENSG00000013016_EHD3_NT.msa.dat", "ENSG00000013016_EHD3_NT.branchlength.dat", "ENSG00000013016_EHD3_NT.table.dat", "5.209", "2", "-9147.9319"),
    ("ENSG00000013016_EHD3_NT.msa.dat", "ENSG00000013016_EHD3_NT.branchlength.dat", "ENSG00000013016_EHD3_NT.table.dat", "5.209", "0.4242", "-6737.3043"),
]

    
## Create a list to store the results
results = []

# Loop through each dataset
for dataset in datasets:
    msa_path, branch_lengths_path, tree_path, branch_len, mu, expected_output = dataset
    print(f"Testing for dataset: {msa_path} (mu={mu})")
    
    # Construct the command to run the code
    command = [
        "python3",
        "run_lkl_calculation.py",
        "--tree_path", f"data/{tree_path}",
        "--branch_lengths_path", f"data/{branch_lengths_path}",
        "--msa_path", f"data/{msa_path}",
        "--mu", mu,
        "--output_lkl_root_node_only"
    ]

    # Run the command and capture the output
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    output = output.decode().strip()

    # Check if the output matches the expected value
    match_status = "Match" if output == expected_output else "No Match"
    
    # Append the result to the list
    results.append([msa_path, mu, expected_output, output, match_status])

# Print the results in a tabular format
headers = ["Dataset", "Mu", "Expected Output", "Actual Output", "Match Status"]
print(tabulate(results, headers=headers, tablefmt="grid"))
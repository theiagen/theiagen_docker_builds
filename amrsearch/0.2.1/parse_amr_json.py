import json
import csv
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Parse AMR JSON output and generate CSV and PNG files.")
parser.add_argument("input_json", help="Path to the input JSON file.")
parser.add_argument("samplename", help="Sample name for output files.")
args = parser.parse_args()

with open(args.input_json, 'r') as f:
    data = json.load(f)

library_version = data['library']['version']
library_label = data['library']['label']
resistance_profiles = data['resistanceProfile']

# Define output filenames using samplename
output_csv = f"{args.samplename}_amr_results.csv"
output_pdf = f"{args.samplename}_amr_results.pdf"
output_version_txt = "output_amr_version.txt"

# Save the AMR search version to output
with open(output_version_txt, 'w') as f:
    f.write(library_version + "\n")

table_data = []
for profile in resistance_profiles:
    agent = profile['agent']['name']
    inferred_resistance = profile['state'].capitalize()  # Capitalizing for consistency
    determinants = []

    acquired = profile['determinants'].get('acquired', [])
    for item in acquired:
        determinants.append(item.get('gene', ''))

    variants = profile['determinants'].get('variants', [])
    for item in variants:
        gene = item.get('gene', '')
        variant = item.get('variant', '')
        determinants.append(f"{gene}_{variant}")

    determinants_str = '; '.join(determinants) if determinants else 'none'
    table_data.append([agent, inferred_resistance, determinants_str])

# Write to CSV
with open(output_csv, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Agent", "Inferred Resistance", "Known Determinants"])
    csvwriter.writerows(table_data)

# Generate PDF with title
fig, ax = plt.subplots(figsize=(16, len(table_data) * 0.5))
ax.set_title(f"AMR - Antimicrobial Resistance (Library {library_label}, Version {library_version})", fontsize=18, weight='bold')
ax.axis('off')

# Create table
table = plt.table(cellText=table_data, colLabels=["Agent", "Inferred Resistance", "Known Determinants"],
                cellLoc='center', loc='center', bbox=[0, 0, 1, 1])

# Style headers
for col, label in enumerate(["Agent", "Inferred Resistance", "Known Determinants"]):
    cell = table[0, col]
    cell.set_text_props(weight="bold", fontsize=17)
    cell.set_facecolor("lightgrey")

# Style rows and highlight "Resistant" rows
for row_idx, row in enumerate(table_data, start=1):
    for col_idx, value in enumerate(row):
        cell = table[row_idx, col_idx]
        cell.set_text_props(weight="normal", fontsize=16)
        if row[1] == "Resistant":  # Highlight row if "Resistant"
            cell.set_text_props(weight="bold", fontsize=16)

plt.savefig(output_pdf, bbox_inches='tight', dpi=300)
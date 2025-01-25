import matplotlib.pyplot as plt
import pandas as pd

# Step 1: Load the data
file_path = "/Users/pavankumar/Downloads/high_impact_genes.txt"  # Update the path if needed
data = pd.read_csv(file_path, sep="\s+", names=["Count", "Gene"])

# Step 2: Sort the data by count (descending order)
sorted_data = data.sort_values(by="Count", ascending=False)

# Step 3: Filter the top N genes for visualization (e.g., top 20)
top_genes = sorted_data.head(20)

# Step 4: Create the bar chart
plt.figure(figsize=(12, 8))
plt.barh(top_genes["Gene"], top_genes["Count"], color="skyblue")
plt.xlabel("Number of High-Impact Variants", fontsize=12)
plt.ylabel("Gene", fontsize=12)
plt.title("Top Genes with High-Impact Variants", fontsize=14)
plt.gca().invert_yaxis()  # Reverse the order to show the highest on top
plt.tight_layout()

# Step 5: Show or save the plot
plt.savefig("high_impact_genes_chart.png")  # Save the plot
plt.show()  # Display the plot

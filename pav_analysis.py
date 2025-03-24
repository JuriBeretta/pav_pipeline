# Define the names for your lists
samples=[]
filepath="/mnt/Archive1/juri/pav_littorina/reads/sample_list.txt"

with open (filepath) as f:
    for line in f:
        samples.append(line.strip())


list_absent_genes = [f"{sample}.sorted.mapping.depth.absent.list" for sample in samples]


ecotypes = ["wave", "crab", "barnacle", "brakish"]


for ecotype in ecotypes:
    genes = set()
    for name in list_absent_genes:
        for sample in samples:
            filepath = f"{ecotype}/{sample}.sorted.mapping.depth.absent.list"
            try:
                with open(filepath) as f:
                    for line in f:
                        genes.add(line.strip())  # Remove duplicates
            except FileNotFoundError:
                print(f"File not found: {filepath}")
            continue
        
    with open(f"merged_{ecotype}.txt", "w") as out:
        for gene in sorted(genes):
                if not gene.startswith ("Trna"): #exclude Trna genes
                    out.write(gene + "\n")

list_merged_files=[f"merged_{ecotype}.txt" for ecotype in ecotypes]

# Load all merged gene lists

all_sets = []
for name in list_merged_files:
    filepath = f"{name}" #A path when needed
    with open(filepath) as f:
        all_sets.append(set(f.read().splitlines()))


# Find unique genes for each list
for i, name in enumerate(list_merged_files):

    # Create a copy of all sets except the current one
    other_sets = all_sets.copy()
    current_set = other_sets.pop(i)
    
    # Union of all other sets
    
    gene_unique = set.union(*other_sets) if other_sets else set()
    

    unique_to_current = current_set - gene_unique
    
    # Write the unique genes to a file
    for ecotype in ecotypes:
        with open(f"unique_to_{ecotype}.txt", "w") as out:
            for gene in sorted(unique_to_current):
                out.write(gene + "\n")
# Find unique genes for each list
for i, name in enumerate(list_merged_files):

    # Create a copy of all sets except the current one
    other_sets = all_sets.copy()
    current_set = other_sets.pop(i)
    
    # Union of all other sets
    
    gene_unique = set.union(*other_sets) if other_sets else set()
    

    unique_to_current = current_set - gene_unique
    
    # Write the unique genes to a file
    for ecotype in ecotypes:
        with open(f"unique_to_{ecotype}.txt", "w") as out:
            for gene in sorted(unique_to_current):
                out.write(gene + "\n")
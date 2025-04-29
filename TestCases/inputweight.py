import random

input_file = "input5.config"  # Your original file
output_file = "input8.config"  # New file with weights

# Range for random weights (adjust as needed)
min_weight = 1
max_weight = 100

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    # Read number of vertices and edges
    V = f_in.readline().strip()
    E = f_in.readline().strip()
    
    # Write V and E to output file
    f_out.write(f"{V}\n")
    f_out.write(f"{E}\n")
    
    # Process each edge and add a random weight
    for line in f_in:
        src, dest, weight = map(int, line.strip().split())
        weight = random.randint(min_weight, max_weight)  # Generate random weight
        f_out.write(f"{src} {dest} {weight}\n")

print(f"Weighted graph written to {output_file}")
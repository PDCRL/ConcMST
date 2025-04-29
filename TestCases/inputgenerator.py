import random
import sys

random.seed(42)  # Optional: for reproducibility

if len(sys.argv) != 2:
    print("Usage: python inputgenerator.py <number_of_vertices>")
    sys.exit(1)

V = int(sys.argv[1])  # number of vertices
E = V * 20            # number of edges

outputFileName = "input" + str(V) + ".config"

edges = set()
edge_pairs = set()

# Step 1: Build a spanning tree (ensures connectivity)
vertices = list(range(V))
random.shuffle(vertices)

for i in range(1, V):
    u = vertices[i]
    v = vertices[random.randint(0, i - 1)]  # connect to an already added vertex
    u_, v_ = min(u, v), max(u, v)
    wt = random.randint(1, 1000)
    edges.add((u_, v_, wt))
    edge_pairs.add((u_, v_))

# Step 2: Add remaining edges randomly
while len(edges) < E:
    u = random.randint(0, V - 1)
    v = random.randint(0, V - 1)
    if u != v:
        u_, v_ = min(u, v), max(u, v)
        if (u_, v_) not in edge_pairs:
            wt = random.randint(1, 1000)
            edges.add((u_, v_, wt))
            edge_pairs.add((u_, v_))

# Step 3: Write to file
with open(outputFileName, "w") as f:
    f.write(f"{V}\n")
    f.write(f"{E}\n")
    sorted_edges = sorted(edges, key=lambda x: (x[0], x[1]))
    for u, v, wt in sorted_edges:
        f.write(f"{u} {v} {wt}\n")
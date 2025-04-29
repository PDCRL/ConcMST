import random

V = 1000  # number of vertices
E = 20000  # number of edges

edges = set()
edge_pairs = set()  # for fast (u, v) uniqueness check

while len(edges) < E:
    u = random.randint(0, V - 1)
    v = random.randint(0, V - 1)
    if u != v:
        u_, v_ = min(u, v), max(u, v)
        if (u_, v_) not in edge_pairs:
            wt = random.randint(1, 1000)
            edges.add((u_, v_, wt))
            edge_pairs.add((u_, v_))

with open("input8.config", "w") as f:
    f.write(f"{V}\n")
    f.write(f"{E}\n")
    sorted_edges = sorted(edges, key=lambda x: (x[0], x[1]))
    for u, v, w in sorted_edges:
        f.write(f"{u} {v} {w}\n")

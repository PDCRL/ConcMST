#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Edge {
public:
    int src, dest, weight;
    Edge() : src(0), dest(0), weight(0) {}
    Edge(int s, int d, int w) : src(s), dest(d), weight(w) {}

    bool operator<(const Edge &e) const {
        return weight < e.weight;
    }
};

class DisjointSet {
    vector<int> parent, rank;

public:
    DisjointSet(int n) {
        parent.resize(n, 0);
        rank.resize(n, 0);
        for (int i = 0; i < n; i++)
            parent[i] = i;
    }

    int find(int u);
    void merge(int x, int y);
};

int DisjointSet::find(int u) {
    if (parent[u] != u)
        parent[u] = find(parent[u]);
    return parent[u];
}

void DisjointSet::merge(int x, int y) {
    int rootX = find(x);
    int rootY = find(y);

    if (rootX != rootY) {
        if (rank[rootX] < rank[rootY])
            parent[rootX] = rootY;
        else if (rank[rootX] > rank[rootY])
            parent[rootY] = rootX;
        else {
            parent[rootY] = rootX;
            rank[rootX]++;
        }
    }
}

class MSTGraph {
    int V, E;
    vector<Edge> edges;

public:
    int mstWeight;
    vector<Edge> mstEdges;

    MSTGraph(int V) : V(V), mstWeight(0) {}
    MSTGraph(int V, int E) : V(V), E(E), mstWeight(0) {}

    void addEdge(int src, int dest, int weight);
    void kruskalMST();
};

void MSTGraph::addEdge(int src, int dest, int weight) {
    edges.push_back(Edge(src, dest, weight));
}

void MSTGraph::kruskalMST() {
    DisjointSet ds(V);

    sort(edges.begin(), edges.end());

    for (const Edge &edge : edges) {
        int set1 = ds.find(edge.src);
        int set2 = ds.find(edge.dest);

        if (set1 != set2) {
            mstEdges.push_back(edge);
            mstWeight += edge.weight;
            ds.merge(set1, set2);

            if (mstEdges.size() == V - 1)
                break;
        }
    }
}

int main() {
    int fileNum = 1;
    cout << "Enter the input file number: ";
    cin >> fileNum;

    string fileName = "input" + to_string(fileNum) + ".config";
    ifstream configFile(fileName);

    if (!configFile) {
        cerr << "Error opening file: " << fileName << endl;
        return 1;
    }

    int k, V, E;
    configFile >> k; // Read the number of threads
    configFile >> V; // Read the number of vertices
    configFile >> E; // Read the number of edges

    MSTGraph kg(V);

    for (int i = 0; i < E; i++) {
        int src, dest, weight;
        configFile >> src >> dest >> weight;
        kg.addEdge(src, dest, weight);
    }

    kg.kruskalMST();
    vector<Edge> mst = kg.mstEdges;

    cout << "Edges in the Minimum Spanning Tree (Kruskal's Algorithm):\n";
    for (const Edge &edge : mst) {
        cout << edge.src << " -- " << edge.dest << " : " << edge.weight << "\n";
    }
    cout << "Total weight of MST: " << kg.mstWeight << endl;

    return 0;
}

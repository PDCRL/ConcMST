#include <chrono>
#include <climits>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Edge {
public:
    int src, dest, weight;
    Edge() : src(-1), dest(-1), weight(numeric_limits<int>::max()) {}
    Edge(int s, int d, int w) : src(s), dest(d), weight(w) {}
};

class DisjointSet {
    vector<int> parent, rank;

public:
    DisjointSet(int n) {
        parent.resize(n);
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
    long long mstWeight;
    vector<Edge> mstEdges;

    MSTGraph(int V) : V(V), mstWeight(0) {}
    MSTGraph(int V, int E) : V(V), E(E), mstWeight(0) {}

    void addEdge(int src, int dest, int weight);
    void boruvkaMST(int numThreads);
};

void MSTGraph::addEdge(int src, int dest, int weight) {
    edges.push_back(Edge(src, dest, weight));
    E++;
}

void MSTGraph::boruvkaMST(int numThreads) {
    int numTrees = V;
    DisjointSet ds(V);
    vector<Edge> cheapest(V, Edge(-1, -1, INT_MAX));

    while (numTrees > 1) {
        fill(cheapest.begin(), cheapest.end(), Edge(-1, -1, INT_MAX));

        for (const Edge &edge : edges) {
            int set1 = ds.find(edge.src);
            int set2 = ds.find(edge.dest);

            if (set1 != set2) {
                if (edge.weight < cheapest[set1].weight)
                    cheapest[set1] = edge;

                if (edge.weight < cheapest[set2].weight)
                    cheapest[set2] = edge;
            }
        }

        vector<Edge> cheapestEdges;
        for (int i = 0; i < V; i++) {
            if (cheapest[i].src != -1) {
                cheapestEdges.push_back(cheapest[i]);
            }
        }

        // Step 3: Process cheapest edges to build MST
        for (const Edge &edge : cheapestEdges) {
            if (edge.src != -1) {
                int set1 = ds.find(edge.src);
                int set2 = ds.find(edge.dest);

                if (set1 != set2) {
                    mstEdges.push_back(edge);
                    mstWeight += edge.weight;
                    ds.merge(set1, set2);
                    numTrees--;
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <file_number>" << " <num_threads>" << endl;
        return 1;
    }

    int fileNum = stoi(argv[1]);
    int k = stoi(argv[2]);
    // cout << "Enter the input file number: ";
    // cin >> fileNum;

    string fileName = "input" + to_string(fileNum) + ".config";
    ifstream configFile(fileName);

    if (!configFile) {
        cerr << "Error opening file: " << fileName << endl;
        return 1;
    }

    int V, E;
    configFile >> V; // Read the number of vertices
    configFile >> E; // Read the number of edges

    MSTGraph bg(V, E);

    for (int i = 0; i < E; i++) {
        int src, dest, weight;
        configFile >> src >> dest >> weight;
        bg.addEdge(src, dest, weight);
    }
    configFile.close();

    auto start = chrono::high_resolution_clock::now();
    bg.boruvkaMST(k);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
    float duration_ms = duration.count() / 1e6f;

    vector<Edge> mst = bg.mstEdges;

    cout << "Edges in the Minimum Spanning Tree (Boruvka's Algorithm):\n";
    for (const Edge &edge : mst) {
        cout << edge.src << " -- " << edge.dest << " : " << edge.weight << "\n";
    }
    cout << "Total weight of MST: " << bg.mstWeight << endl;
    cout << "Time taken to build MST: " << duration.count() << " nanoseconds" << endl;
    cout << "Time taken to build MST: " << duration_ms << " milliseconds" << endl;

    return 0;
}

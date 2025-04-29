#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Edge {
public:
    int src, dest, weight;
    Edge() : src(0), dest(0), weight(0) {}
    Edge(int s, int d, int w) : src(s), dest(d), weight(w) {}
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
    long long mstWeight;
    vector<Edge> mstEdges;

    MSTGraph(int V) : V(V), E(0), mstWeight(0) {}
    MSTGraph(int V, int E) : V(V), E(E), mstWeight(0) {}

    void addEdge(int src, int dest, int weight);
    void boruvkaMST(int numThreads);
};

void MSTGraph::addEdge(int src, int dest, int weight) {
    edges.push_back(Edge(src, dest, weight));
}

void MSTGraph::boruvkaMST(int numThreads) {
    DisjointSet ds(V);
    vector<Edge> cheapest(V);
    int numTrees = V;

    while (numTrees > 1) {
        for (int i = 0; i < V; i++)
            cheapest[i] = Edge(-1, -1, INT_MAX);

        for (const Edge &edge : edges) {
            int set1 = ds.find(edge.src);
            int set2 = ds.find(edge.dest);
            cout << "Processing edge: " << edge.src << " -- " << edge.dest << " : " << edge.weight << endl;
            cout << "Set1: " << set1 << ", Set2: " << set2 << endl;

            if (set1 != set2) {
                if (edge.weight < cheapest[set1].weight) {
                    cout << "Cheapest edge for set " << set1 << " updated from " << cheapest[set1].src << " -- " << cheapest[set1].dest << " : " << cheapest[set1].weight << " to: " << edge.src << " -- " << edge.dest << " : " << edge.weight << endl;
                    cheapest[set1] = Edge(edge.src, edge.dest, edge.weight);
                }

                if (edge.weight < cheapest[set2].weight) {
                    cout << "Cheapest edge for set " << set2 << " updated from " << cheapest[set2].src << " -- " << cheapest[set2].dest << " : " << cheapest[set2].weight << " to: " << edge.src << " -- " << edge.dest << " : " << edge.weight << endl;
                    cheapest[set2] = Edge(edge.src, edge.dest, edge.weight);
                }
            }
            cout << "---------------------------------------" << endl;
        }

        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        for (int i = 0; i < V; i++) {
            cout << "Cheapest edge for set " << i << ": " << cheapest[i].src << " -- " << cheapest[i].dest << " : " << cheapest[i].weight << endl;
        }
        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;

        for (int i = 0; i < V; i++) {
            Edge edge = cheapest[i];
            if (edge.src != -1) {
                int set1 = ds.find(edge.src);
                int set2 = ds.find(edge.dest);

                cout << "Processing cheapest edge: " << edge.src << " -- " << edge.dest << " : " << edge.weight << endl;
                cout << "Set1: " << set1 << ", Set2: " << set2 << endl;

                if (set1 != set2) {
                    mstEdges.push_back(edge);
                    mstWeight += edge.weight;
                    ds.merge(set1, set2);

                    cout << "Edge added to MST: " << edge.src << " -- " << edge.dest << " : " << edge.weight << endl;

                    numTrees--;
                    cout << "Number of trees: " << numTrees << endl;
                    cout << "#######################################" << endl;
                }
            }
        }
        cout << "=======================================" << endl;
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

    MSTGraph bg(V, E);

    for (int i = 0; i < E; i++) {
        int src, dest, weight;
        configFile >> src >> dest >> weight;
        bg.addEdge(src, dest, weight);
    }

    bg.boruvkaMST(k);
    vector<Edge> mst = bg.mstEdges;

    cout << "Edges in the Minimum Spanning Tree (Boruvka's Algorithm):\n";
    for (const Edge &edge : mst) {
        cout << edge.src << " -- " << edge.dest << " : " << edge.weight << "\n";
    }
    cout << "Total weight of MST: " << bg.mstWeight << endl;

    return 0;
}

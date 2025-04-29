#include <climits>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

using namespace std;

class Edge {
public:
    int src, dest, weight;
    Edge() : src(0), dest(0), weight(0) {}
    Edge(int s, int d, int w) : src(s), dest(d), weight(w) {}
};

class MSTGraph {
    int V, E;
    vector<vector<pair<int, int>>> adj;

public:
    int mstWeight;
    vector<Edge> mstEdges;

    MSTGraph(int V) : V(V), mstWeight(0) {
        adj.resize(V);
    }

    MSTGraph(int V, int E) : V(V), E(E), mstWeight(0) {
        adj.resize(V);
    }

    void addEdge(int src, int dest, int weight);
    void primMST();
};

void MSTGraph::addEdge(int src, int dest, int weight) {
    adj[src].push_back({dest, weight});
    adj[dest].push_back({src, weight});
}

void MSTGraph::primMST() {
    vector<bool> inMST(V, false);
    vector<int> key(V, INT_MAX);

    auto cmp = [](Edge a, Edge b) { return a.weight > b.weight; };
    priority_queue<Edge, vector<Edge>, decltype(cmp)> pq(cmp);

    key[0] = 0;
    pq.push(Edge(-1, 0, 0));
    mstWeight = 0;

    while (!pq.empty()) {
        Edge e = pq.top();
        auto [u, v, weight] = tuple{e.src, e.dest, e.weight};

        pq.pop();

        if (inMST[v])
            continue; // Skip if vertex v is already in MST as it will form a cycle

        inMST[v] = true;

        if (u != -1) {
            mstEdges.push_back(Edge(u, v, weight));
            mstWeight += weight;
        }

        for (auto &neighbor : adj[v]) {
            int p = v;
            auto [q, w] = neighbor;

            if (!inMST[q] && w < key[q]) {
                key[q] = w;
                pq.push(Edge(p, q, w));
            }
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

    MSTGraph pg(V, E);

    for (int i = 0; i < E; i++) {
        int src, dest, weight;
        configFile >> src >> dest >> weight;
        pg.addEdge(src, dest, weight);
    }

    pg.primMST();
    vector<Edge> mst = pg.mstEdges;

    cout << "Edges in the Minimum Spanning Tree (Prim's Algorithm):\n";
    for (const Edge &edge : mst) {
        cout << edge.src << " -- " << edge.dest << " : " << edge.weight << "\n";
    }
    cout << "Total weight of MST: " << pg.mstWeight << endl;

    return 0;
}

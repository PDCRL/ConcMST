#include <algorithm>
#include <chrono>
#include <climits>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
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
    bool merge(int x, int y);
};

int DisjointSet::find(int u) {
    if (parent[u] != u)
        parent[u] = find(parent[u]);
    return parent[u];
}

bool DisjointSet::merge(int x, int y) {
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
        return true;
    }
    return false;
}

class MSTGraph {
    int V, E;
    vector<Edge> edges;
    vector<mutex> componentLocks; // Per-component locks for cheapest edges
    mutex mstMutex;               // For MST updates

public:
    long long mstWeight;
    vector<Edge> mstEdges;

    MSTGraph(int V) : V(V), E(0), mstWeight(0), componentLocks(V) {}
    MSTGraph(int V, int E) : V(V), E(E), mstWeight(0), componentLocks(V) {}

    void addEdge(int src, int dest, int weight);
    void boruvkaMST(int numThreads);

private:
    void findCheapestEdges(int thread_id, size_t start, size_t end, const vector<Edge> &edges, DisjointSet &ds, vector<Edge> &cheapest);
    void processCheapestEdges(int thread_id, size_t start, size_t end, const vector<Edge> &cheapestEdges, DisjointSet &ds, int &numTrees);
};

void MSTGraph::addEdge(int src, int dest, int weight) {
    edges.push_back(Edge(src, dest, weight));
    E++;
}

void MSTGraph::findCheapestEdges(int thread_id, size_t start, size_t end,
                                 const vector<Edge> &edges, DisjointSet &ds,
                                 vector<Edge> &cheapest) {
    vector<Edge> local_cheapest(V, Edge(-1, -1, INT_MAX));

    for (size_t i = start; i < end; i++) {
        const Edge &edge = edges[i];
        int set1 = ds.find(edge.src);
        int set2 = ds.find(edge.dest);
        if (set1 != set2) {
            if (edge.weight < local_cheapest[set1].weight)
                local_cheapest[set1] = edge;
            if (edge.weight < local_cheapest[set2].weight)
                local_cheapest[set2] = edge;
        }
    }

    for (int v = 0; v < V; v++) {
        if (local_cheapest[v].weight != INT_MAX) {
            lock_guard<mutex> lock(componentLocks[v]);
            if (local_cheapest[v].weight < cheapest[v].weight)
                cheapest[v] = local_cheapest[v];
        }
    }
}

void MSTGraph::processCheapestEdges(int thread_id, size_t start, size_t end, const vector<Edge> &cheapestEdges, DisjointSet &ds, int &numTrees) {
    for (size_t i = start; i < end && i < cheapestEdges.size(); i++) {
        const Edge &edge = cheapestEdges[i];
        if (edge.src != -1) {
            int set1 = ds.find(edge.src);
            int set2 = ds.find(edge.dest);
            if (set1 != set2) {
                lock_guard<mutex> lock(mstMutex);
                if (ds.merge(set1, set2)) {
                    mstEdges.push_back(edge);
                    mstWeight += edge.weight;
                    numTrees--;
                }
            }
        }
    }
}

void MSTGraph::boruvkaMST(int numThreads) {
    int numTrees = V;
    DisjointSet ds(V);
    vector<Edge> cheapest(V, Edge(-1, -1, INT_MAX));

    while (numTrees > 1) {
        fill(cheapest.begin(), cheapest.end(), Edge(-1, -1, INT_MAX));

        vector<thread> threads;
        size_t edgeCount = edges.size();
        size_t chunkSize = (edgeCount + numThreads - 1) / numThreads;

        for (int t = 0; t < numThreads; t++) {
            size_t start = t * chunkSize;
            size_t end = (t == numThreads - 1) ? edgeCount : start + chunkSize;
            threads.emplace_back(&MSTGraph::findCheapestEdges, this, t, start, end, ref(edges), ref(ds), ref(cheapest));
        }

        for (auto &t : threads) {
            t.join();
        }

        vector<Edge> cheapestEdges;
        for (int i = 0; i < V; i++) {
            if (cheapest[i].src != -1) {
                cheapestEdges.push_back(cheapest[i]);
            }
        }

        // Step 3: Parallel cheapest edge processing
        threads.clear();
        size_t numCheapestEdges = cheapestEdges.size();
        chunkSize = (numCheapestEdges + numThreads - 1) / numThreads;

        for (int t = 0; t < numThreads; t++) {
            size_t start = t * chunkSize;
            size_t end = (t == numThreads - 1) ? numCheapestEdges : start + chunkSize;
            threads.emplace_back(&MSTGraph::processCheapestEdges, this, t, start, end, ref(cheapestEdges), ref(ds), ref(numTrees));
        }

        for (auto &t : threads) {
            t.join();
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;

class Edge {
public:
    int src, dest, weight;
    Edge() : src(-1), dest(-1), weight(numeric_limits<int>::max()) {}
    Edge(int s, int d, int w) : src(s), dest(d), weight(w) {}

    // Operator for sorting edges: by src, then by dest itself.
    bool operator<(const Edge &other) const {
        if (src != other.src)
            return src < other.src;
        return dest < other.dest;
    }
};

class Graph {
    int V, E;
    unordered_map<int, vector<Edge>> adj;

public:
    Graph(int V) : V(V), E(0) {}
    Graph(int V, int E) : V(V), E(E) {}

    void addEdge(int src, int dest, int weight);
    vector<Edge> getEdgeList();
};

void Graph::addEdge(int src, int dest, int weight) {
    adj[src].push_back(Edge(src, dest, weight));
    adj[dest].push_back(Edge(dest, src, weight));
    E++;
}

vector<Edge> Graph::getEdgeList() {
    vector<Edge> edgeList;

    for (const auto &vertex : adj) {
        int src = vertex.first;
        for (const Edge &edge : vertex.second) {
            if (src < edge.dest) { // Process edge only once
                edgeList.push_back(edge);
            }
        }
    }

    return edgeList;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <file_number>" << endl;
        return 1;
    }

    int fileNum = stoi(argv[1]);
    string inputFileName = "input" + to_string(fileNum) + ".config";
    ifstream configFile(inputFileName);

    if (!configFile) {
        cerr << "Error opening input file: " << inputFileName << endl;
        return 1;
    }

    int V, E;
    configFile >> V; // Read the number of vertices
    configFile >> E; // Read the number of edges

    Graph g(V, E);

    for (int i = 0; i < E; i++) {
        int src, dest, weight;
        configFile >> src >> dest >> weight;
        g.addEdge(src, dest, weight);
    }
    configFile.close();

    auto start = chrono::high_resolution_clock::now();
    vector<Edge> edgeList = g.getEdgeList();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

    // Open output file
    string outputFileName = "output" + to_string(fileNum) + ".txt";
    ofstream outputFile(outputFileName);
    if (!outputFile) {
        cerr << "Error opening output file: " << outputFileName << endl;
        return 1;
    }

    int numEdges = edgeList.size();
    outputFile << fileNum << endl;  // Write file number to output file
    outputFile << numEdges << endl; // Write number of edges to output file

    // Sort edges by weight, then by src, then by dest
    sort(edgeList.begin(), edgeList.end());

    for (const Edge &edge : edgeList) {
        outputFile << edge.src << " " << edge.dest << " " << edge.weight << "\n";
    }

    cout << "Time taken to get edge list: " << duration.count() << " microseconds" << endl;
    outputFile.close();
    cout << "Sorted edge list written to: " << outputFileName << endl;

    return 0;
}
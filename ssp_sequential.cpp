#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#include <climits>
#include <chrono>
#include <random>  // For random edge generation
#include <algorithm>  // For shuffle function

using namespace std;
using namespace std::chrono;

const int MAX_DISTANCE = INT_MAX;

struct GraphEdge {
    int to;
    int weight;
};

struct EdgeChange {
    int from;
    int to;
    int weight;
    bool isAddition;
};

vector<vector<GraphEdge>> adjacencyList;
unordered_map<int, int> idToIndexMap;
vector<int> indexToId;
int nextIndex = 0;
int totalNodes = 0;

vector<int> distances, predecessors;
vector<bool> needsUpdate;

// Load graph data from an unweighted edge list file
void loadGraphFromFile(const string& filename) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error opening input file.\n";
        exit(1);
    }

    string line;
    int total_edges = 0;
    while (getline(infile, line)) {
        if (line.size() == 0 || line[0] == '#') continue;  // Skip comments or empty lines
        stringstream ss(line);
        int u_orig, v_orig;
        ss >> u_orig >> v_orig;

        // Assign indices to new nodes
        if (!idToIndexMap.count(u_orig)) {
            idToIndexMap[u_orig] = nextIndex++;
            indexToId.push_back(u_orig);
        }
        if (!idToIndexMap.count(v_orig)) {
            idToIndexMap[v_orig] = nextIndex++;
            indexToId.push_back(v_orig);
        }

        int u = idToIndexMap[u_orig];
        int v = idToIndexMap[v_orig];

        if (adjacencyList.size() <= max(u, v)) {
            adjacencyList.resize(max(u, v) + 1);
        }

        adjacencyList[u].push_back({v, 1});  // weight = 1 for unweighted graph
        total_edges++;
    }

    totalNodes = nextIndex;
    distances.resize(totalNodes, MAX_DISTANCE);
    predecessors.resize(totalNodes, -1);
    needsUpdate.resize(totalNodes, false);
    cout << "Loaded " << totalNodes << " nodes and " << total_edges << " edges from input file.\n";
}

// Prepare data structures for Dijkstra
void initializeSingleSource(int src) {
    fill(distances.begin(), distances.end(), MAX_DISTANCE);
    fill(predecessors.begin(), predecessors.end(), -1);
    distances[src] = 0;
}

// Compute shortest paths with Dijkstra
void runDijkstra(int src) {
    auto startTime = high_resolution_clock::now();

    initializeSingleSource(src);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.push({0, src});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > distances[u]) continue;

        for (const GraphEdge& edge : adjacencyList[u]) {
            int v = edge.to;
            if (distances[v] > distances[u] + edge.weight) {
                distances[v] = distances[u] + edge.weight;
                predecessors[v] = u;
                pq.push({distances[v], v});
            }
        }
    }

    auto endTime = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(endTime - startTime).count();
    cout << "Dijkstra algorithm completed in " << duration << " ms\n";
}

// Check if an edge exists in the graph
bool hasEdge(int u, int v) {
    for (const GraphEdge& edge : adjacencyList[u]) {
        if (edge.to == v) {
            return true;
        }
    }
    return false;
}

// Add a new node to the graph if not present
void ensureNodeExists(int nodeId) {
    if (idToIndexMap.find(nodeId) == idToIndexMap.end()) {
        idToIndexMap[nodeId] = nextIndex++;
        indexToId.push_back(nodeId);
        adjacencyList.resize(nextIndex);
        distances.push_back(MAX_DISTANCE);
        predecessors.push_back(-1);
        needsUpdate.push_back(false);
    }
}

// Apply a single edge update and then recompute shortest paths
bool updateEdgeAndRecompute(const EdgeChange& update, int sourceIndex) {
    auto startTime = high_resolution_clock::now();
    bool success = true;

    if (update.isAddition) {
        if (hasEdge(update.from, update.to)) {
            cout << "Edge already exists: " << indexToId[update.from] << " -> " 
                 << indexToId[update.to] << "\n";
            success = false;
        } else {
            cout << "Inserting edge: " << indexToId[update.from] << " -> " 
                 << indexToId[update.to] << " (weight: " << update.weight << ")\n";
            adjacencyList[update.from].push_back({update.to, update.weight});
        }
    } else {
        if (!hasEdge(update.from, update.to)) {
            cout << "Edge doesn't exist: " << indexToId[update.from] << " -> " 
                 << indexToId[update.to] << "\n";
            success = false;
        } else {
            cout << "Removing edge: " << indexToId[update.from] << " -> " 
                 << indexToId[update.to] << "\n";
            auto target = update.to;
            auto it = find_if(adjacencyList[update.from].begin(), adjacencyList[update.from].end(),
                              [target](const GraphEdge& edge) { return edge.to == target; });
            if (it != adjacencyList[update.from].end()) {
                adjacencyList[update.from].erase(it);
            }
        }
    }

    // Only recompute if the graph was modified
    if (success) {
        // Recompute shortest paths after each update
        runDijkstra(sourceIndex);
        return true;
    }
    return false;
}

// Apply a batch of edge updates (incremental)
int applyBatchUpdates(vector<EdgeChange>& updates) {
    int affectedNodesCount = 0;

    for (auto& update : updates) {
        if (update.isAddition) {
            adjacencyList[update.from].push_back({update.to, update.weight});
            if (distances[update.from] != MAX_DISTANCE && distances[update.to] > distances[update.from] + update.weight) {
                distances[update.to] = distances[update.from] + update.weight;
                predecessors[update.to] = update.from;
                needsUpdate[update.to] = true;
                affectedNodesCount++;
            }
        } else {
            auto target = update.to;
            auto it = find_if(adjacencyList[update.from].begin(), adjacencyList[update.from].end(),
                              [target](const GraphEdge& edge) { return edge.to == target; });
            if (it != adjacencyList[update.from].end()) {
                adjacencyList[update.from].erase(it);
            }
            if (predecessors[update.to] == update.from) {
                distances[update.to] = MAX_DISTANCE;
                predecessors[update.to] = -1;
                needsUpdate[update.to] = true;
                affectedNodesCount++;
            }
        }
    }

    cout << "Number of nodes affected by changes: " << affectedNodesCount << " (" 
         << (float)affectedNodesCount / totalNodes * 100 << "% of total nodes)\n";

    return affectedNodesCount;
}

// Iteratively update distances for nodes flagged for updates
void updateAffectedDistances() {
    auto startTime = high_resolution_clock::now();

    int iterations = 0;
    int totalNodesUpdated = 0;

    bool distancesChanged = true;
    while (distancesChanged) {
        distancesChanged = false;
        iterations++;
        int updatesInIteration = 0;

        for (int u = 0; u < totalNodes; ++u) {
            if (!needsUpdate[u]) continue;
            needsUpdate[u] = false;

            for (const GraphEdge& edge : adjacencyList[u]) {
                int v = edge.to;
                if (distances[v] > distances[u] + edge.weight) {
                    distances[v] = distances[u] + edge.weight;
                    predecessors[v] = u;
                    needsUpdate[v] = true;
                    distancesChanged = true;
                    updatesInIteration++;
                }
            }
        }

        totalNodesUpdated += updatesInIteration;
        cout << "  - Iteration " << iterations << ": Updated " << updatesInIteration << " nodes\n";
    }

    auto endTime = high_resolution_clock::now();
    auto durationMs = duration_cast<milliseconds>(endTime - startTime).count();
    cout << "Incremental update completed in " << durationMs << " ms\n";
    cout << "Total iterations: " << iterations << ", Total node updates: " << totalNodesUpdated << "\n";
}

// Generate random edge updates (insertions or deletions) for the requested count
vector<EdgeChange> generateRandomUpdates(int numUpdates, bool isInsertion) {
    vector<EdgeChange> updates;

    // Random number generator setup
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> nodeDist(0, totalNodes - 1);

    if (isInsertion) {
        // For insertion mode, avoid iterating all possible edges
        // Generate random non-existing edges until required count is reached
        int attemptCount = 0;
        int maxAttempts = numUpdates * 100;
        cout << "Generating " << numUpdates << " random insertions...\n";
        while (updates.size() < numUpdates && attemptCount < maxAttempts) {
            attemptCount++;

            // Pick two random nodes
            int u = nodeDist(gen);
            int v = nodeDist(gen);

            // Skip if it's a self-loop (same node)
            if (u == v) continue;

            // Skip if the edge already exists
            if (!hasEdge(u, v)) {
                updates.push_back({u, v, 1, true});
                if (updates.size() % 10 == 0 || updates.size() == numUpdates) {
                    cout << "  Progress: " << updates.size() << "/" << numUpdates << " edges found\r";
                    cout.flush();
                }
            }
        }
        cout << endl;

        if (updates.size() < numUpdates) {
            cout << "Warning: Could only generate " << updates.size() << " valid insertions after " 
                 << attemptCount << " attempts.\n";
        }
    } else {
        // For deletion mode, select random existing edges
        vector<pair<int, int>> allEdges;

        cout << "Collecting all existing edges...\n";
        // Collect all existing edges
        for (int u = 0; u < totalNodes; u++) {
            for (const auto& edge : adjacencyList[u]) {
                allEdges.push_back({u, edge.to});
            }
        }

        cout << "Found " << allEdges.size() << " existing edges.\n";

        // If not enough edges are available, adjust the count
        if (allEdges.size() < numUpdates) {
            cout << "Warning: only " << allEdges.size() << " edges available for deletion; using all of them.\n";
            numUpdates = allEdges.size();
        }

        // Shuffle the list of existing edges
        cout << "Randomizing edge order...\n";
        shuffle(allEdges.begin(), allEdges.end(), gen);

        // Take the first numUpdates edges for deletion
        for (int i = 0; i < numUpdates; i++) {
            updates.push_back({allEdges[i].first, allEdges[i].second, 1, false});
        }
    }

    cout << "Generated " << updates.size() << " valid " 
         << (isInsertion ? "insertions" : "deletions") << " successfully.\n";

    return updates;
}

// Display shortest paths from source (distance and predecessor)
void printShortestPaths(int limit = 100) {
    cout << "\nShortest paths from source:\n";
    for (int i = 0; i < min(totalNodes, limit); ++i) {
        int distVal = (distances[i] == MAX_DISTANCE ? -1 : distances[i]);
        int predVal = (predecessors[i] == -1 ? -1 : indexToId[predecessors[i]]);
        cout << "Node " << indexToId[i]
             << ": Distance = " << distVal
             << ", Predecessor = " << predVal << endl;
    }
}

int main() {
    loadGraphFromFile("graph.txt");

    int sourceId;
    cout << "Enter source node ID: ";
    cin >> sourceId;

    if (!idToIndexMap.count(sourceId)) {
        cerr << "Invalid source node ID.\n";
        return 1;
    }

    int sourceIndex = idToIndexMap[sourceId];
    runDijkstra(sourceIndex);
    printShortestPaths();

    int menuOption;
    do {
        cout << "\n1. Random Batch Updates\n2. Manual Updates\n3. Exit\nChoice: ";
        cin >> menuOption;

        if (menuOption == 1) {
            // Random batch updates
            int numUpdates;
            int actionType;
            cout << "Enter number of updates: ";
            cin >> numUpdates;

            cout << "Update type (1: Insert, 2: Delete): ";
            cin >> actionType;

            bool isInsertion = (actionType == 1);

            vector<EdgeChange> updates = generateRandomUpdates(numUpdates, isInsertion);

            cout << "\nGenerating " << numUpdates << " random " 
                 << (isInsertion ? "insertions" : "deletions") << "...\n";

            auto totalStart = high_resolution_clock::now();

            // Apply each update one by one, with full recomputation each time
            int successfulUpdates = 0;

            for (const auto& update : updates) {
                auto updateStartTime = high_resolution_clock::now();
                bool success = updateEdgeAndRecompute(update, sourceIndex);
                auto updateEndTime = high_resolution_clock::now();

                if (success) {
                    successfulUpdates++;
                }
            }

            auto totalEnd = high_resolution_clock::now();
            auto totalDuration = duration_cast<milliseconds>(totalEnd - totalStart).count();

            cout << "\nBatch update summary:\n";
            cout << "Successful updates: " << successfulUpdates << " out of " << updates.size() << "\n";
            cout << "Total time: " << totalDuration << " ms\n";
            if (successfulUpdates > 0) {
                cout << "Average time per successful update: " << totalDuration / successfulUpdates << " ms\n";
            }

            printShortestPaths();
        } else if (menuOption == 2) {
            // Manual updates
            int modeChoice;
            cout << "\n1. Insert Batch\n2. Delete Batch\nChoice: ";
            cin >> modeChoice;

            if (modeChoice == 1 || modeChoice == 2) {
                int numUpdates;
                cout << "Enter number of changes: ";
                cin >> numUpdates;
                vector<EdgeChange> updates;

                for (int i = 0; i < numUpdates; ++i) {
                    int fromId, toId;
                    cout << "EdgeChange " << i + 1 << " (from to): ";
                    cin >> fromId >> toId;

                    ensureNodeExists(fromId);
                    ensureNodeExists(toId);
                    updates.push_back({idToIndexMap[fromId], idToIndexMap[toId], 1, modeChoice == 1});
                }

                auto totalStart = high_resolution_clock::now();

                // Apply each update one by one, with full recomputation each time
                int successfulUpdates = 0;

                for (const auto& update : updates) {
                    auto updateStartTime = high_resolution_clock::now();
                    bool success = updateEdgeAndRecompute(update, sourceIndex);
                    auto updateEndTime = high_resolution_clock::now();

                    if (success) {
                        successfulUpdates++;
                    }
                }

                auto totalEnd = high_resolution_clock::now();
                auto totalDuration = duration_cast<milliseconds>(totalEnd - totalStart).count();

                cout << "\nBatch update summary:\n";
                cout << "Successful updates: " << successfulUpdates << " out of " << updates.size() << "\n";
                cout << "Total time: " << totalDuration << " ms\n";
                if (successfulUpdates > 0) {
                    cout << "Average time per successful update: " << totalDuration / successfulUpdates << " ms\n";
                }

                printShortestPaths();
            }
        }
    } while (menuOption != 3);

    return 0;
}
#include <iostream>
#include <fstream>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <mpi.h>
#include <metis.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace std::chrono;

const int INF = INT_MAX;

struct GraphNode {
    int dest;
    int weight;
    GraphNode* next;
};

struct SSSPGraph {
    int V;
    vector<GraphNode*> adjacencyList;

    SSSPGraph(int vertices) : V(vertices), adjacencyList(vertices, nullptr) {}

    void addEdge(int src, int dest, int weight) {
        removeEdge(src, dest);
        adjacencyList[src] = new GraphNode{dest, weight, adjacencyList[src]};
    }

    void removeEdge(int src, int dest) {
        GraphNode* curr = adjacencyList[src];
        GraphNode* prev = nullptr;
        while (curr) {
            if (curr->dest == dest) {
                if (prev) prev->next = curr->next;
                else adjacencyList[src] = curr->next;
                delete curr;
                break;
            }
            prev = curr;
            curr = curr->next;
        }
    }

    void clear() {
        for (auto& node : adjacencyList) {
            while (node) {
                GraphNode* temp = node;
                node = node->next;
                delete temp;
            }
        }
    }
};

void syncBoundaryNodesAcrossProcesses(SSSPGraph& g, vector<int>& distanceFromSource, vector<int>& predecessor, 
                         const vector<int>& partitionMap, int rank, int size) {
    vector<int> boundary_nodes;
    for (int u = 0; u < g.V; u++) {
        if (partitionMap[u] != rank) continue;
        
        bool is_boundary = false;
        for (GraphNode* curr = g.adjacencyList[u]; curr; curr = curr->next) {
            if (partitionMap[curr->dest] != rank) {
                is_boundary = true;
                break;
            }
        }
        
        if (is_boundary) {
            boundary_nodes.push_back(u);
        }
    }

    int local_count = boundary_nodes.size();
    vector<int> boundary_counts(size);
    
    for (int i = 0; i < size; i++) {
        int count_to_share = (i == rank) ? local_count : 0;
        MPI_Bcast(&count_to_share, 1, MPI_INT, i, MPI_COMM_WORLD);
        boundary_counts[i] = count_to_share;
    }
    
    for (int i = 0; i < size; i++) {
        if (boundary_counts[i] > 0) {
            vector<int> proc_data;
            
            if (i == rank) {
                for (int u : boundary_nodes) {
                    proc_data.push_back(u);
                    proc_data.push_back(distanceFromSource[u]);
                    proc_data.push_back(predecessor[u]);
                }
            } else {
                proc_data.resize(boundary_counts[i] * 3);
            }
            
            MPI_Bcast(proc_data.data(), proc_data.size(), MPI_INT, i, MPI_COMM_WORLD);
            
            if (i != rank) {
                for (size_t j = 0; j < proc_data.size(); j += 3) {
                    int u = proc_data[j];
                    int new_distanceFromSource = proc_data[j+1];
                    int new_predecessor = proc_data[j+2];
                    
                    if (new_distanceFromSource < distanceFromSource[u]) {
                        distanceFromSource[u] = new_distanceFromSource;
                        predecessor[u] = new_predecessor;
                    }
                }
            }
        }
    }
}

void runDistributedDijkstra(SSSPGraph& g, vector<int>& distanceFromSource, vector<int>& predecessor, 
                           const vector<int>& partitionMap, int rank, int size, 
                           long long& total_sssp_time) {
    int updated = 1;
    int iterations = 0;
    const int MAX_ITERATIONS = 1000;
    
    auto start_time = high_resolution_clock::now();
    
    if (rank == 0) {
        cout << "Starting Dijkstra's algorithm..." << endl;
    }
    
    while (updated && iterations++ < MAX_ITERATIONS) {
        updated = 0;
        
        for (int u = 0; u < g.V; u++) {
            if (partitionMap[u] != rank || distanceFromSource[u] == INF) continue;
            
            for (GraphNode* curr = g.adjacencyList[u]; curr; curr = curr->next) {
                int v = curr->dest;
                int new_distanceFromSource = distanceFromSource[u] + curr->weight;
                
                if (new_distanceFromSource < distanceFromSource[v]) {
                    distanceFromSource[v] = new_distanceFromSource;
                    predecessor[v] = u;
                    updated = 1;
                }
            }
        }
        
        syncBoundaryNodesAcrossProcesses(g, distanceFromSource, predecessor, partitionMap, rank, size);
        
        if (rank == 0 && iterations % 10 == 0) {
            cout << "Iteration " << iterations << " completed." << endl;
        }
        
        int global_updated = 0;
        int local_updated = updated;
        MPI_Allreduce(&local_updated, &global_updated, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        updated = global_updated;
    }
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    total_sssp_time += duration.count();
    
    if (rank == 0) {
        cout << "Dijkstra's algorithm completed in " << iterations << " iterations." << endl;
        cout << "This execution took: " << duration.count() << " milliseconds" << endl;
    }
    
    if (iterations >= MAX_ITERATIONS && rank == 0) {
        cerr << "Warning: SSSP computation reached maximum iterations" << endl;
    }
}

void writeProcessOutputToFile(const vector<int>& distanceFromSource, const vector<int>& predecessor, 
                 const vector<int>& indexToOriginalId, int src, int rank) {
    ofstream fout("output_" + to_string(rank) + ".txt");
    fout << "Process " << rank << " results:\n";
    
    int reachable_nodes = 0;
    int total_distanceFromSourceance = 0;
    int max_distanceFromSourceance = 0;
    
    for (size_t i = 0; i < distanceFromSource.size(); i++) {
        if (distanceFromSource[i] != INF) {
            reachable_nodes++;
            total_distanceFromSourceance += distanceFromSource[i];
            max_distanceFromSourceance = max(max_distanceFromSourceance, distanceFromSource[i]);
            
            fout << "Node " << indexToOriginalId[i] << ": Distance = " << distanceFromSource[i] << ", Path = ";
            
            vector<int> path;
            for (int v = i; v != -1; v = predecessor[v]) {
                path.push_back(indexToOriginalId[v]);
            }
            
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                fout << *it << " ";
            }
            fout << "\n";
        }
    }
    
    fout << "\nSummary Statistics:\n";
    fout << "Reachable nodes: " << reachable_nodes << " out of " << distanceFromSource.size() << "\n";
    if (reachable_nodes > 0) {
        fout << "Average distanceFromSourceance: " << fixed << setprecision(2) << (double)total_distanceFromSourceance / reachable_nodes << "\n";
        fout << "Maximum distanceFromSourceance: " << max_distanceFromSourceance << "\n";
    }
    
    fout.close();
    
    int global_reachable = 0;
    int global_total_distanceFromSource = 0;
    int global_max_distanceFromSource = 0;
    
    MPI_Reduce(&reachable_nodes, &global_reachable, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_distanceFromSourceance, &global_total_distanceFromSource, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&max_distanceFromSourceance, &global_max_distanceFromSource, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        cout << "\n===== Dijkstra's Algorithm Results =====" << endl;
        cout << "Source node: " << indexToOriginalId[src] << endl;
        cout << "Reachable nodes: " << global_reachable << " out of " << distanceFromSource.size() << endl;
        if (global_reachable > 0) {
            cout << "Average distanceFromSourceance: " << fixed << setprecision(2) 
                 << (double)global_total_distanceFromSource / global_reachable << endl;
            cout << "Maximum distanceFromSourceance: " << global_max_distanceFromSource << endl;
        }
        cout << "Detailed results written to output_[rank].txt files" << endl;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto overall_start_time = high_resolution_clock::now();
    long long total_sssp_time = 0;
    long long update_time = 0;
    long long graph_read_time = 0;
    long long partitionMapition_time = 0;

    map<int, int> originalIdToIndex;
    vector<int> indexToOriginalId;
    SSSPGraph g(0);

    // SSSPGraph reading
    auto graph_read_start = high_resolution_clock::now();
    if (rank == 0) {
        cout << "Reading graph file..." << endl;
        ifstream fin("graph.txt");
        if (!fin) {
            cerr << "Error: Could not open graph.txt" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        string line;
        int from, to;
        
        while (getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            sscanf(line.c_str(), "%d%d", &from, &to);
            if (originalIdToIndex.find(from) == originalIdToIndex.end()) {
                originalIdToIndex[from] = indexToOriginalId.size();
                indexToOriginalId.push_back(from);
            }
            if (originalIdToIndex.find(to) == originalIdToIndex.end()) {
                originalIdToIndex[to] = indexToOriginalId.size();
                indexToOriginalId.push_back(to);
            }
        }
        
        g = SSSPGraph(indexToOriginalId.size());
        fin.clear();
        fin.seekg(0);
        
        while (getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            int from, to;
            sscanf(line.c_str(), "%d%d", &from, &to);
            g.addEdge(originalIdToIndex[from], originalIdToIndex[to], 1);
        }
        
        auto graph_read_end = high_resolution_clock::now();
        graph_read_time = duration_cast<milliseconds>(graph_read_end - graph_read_start).count();
        cout << "SSSPGraph reading completed. Time taken: " << graph_read_time << " ms" << endl;
        cout << "Number of vertices: " << g.V << endl;
    }

    int num_vertices = 0;
    if (rank == 0) {
        num_vertices = g.V;
    }
    MPI_Bcast(&num_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        g = SSSPGraph(num_vertices);
        indexToOriginalId.resize(num_vertices);
    }
    
    if (num_vertices > 0) {
        vector<int> temp_id_buffer;
        if (rank == 0) {
            temp_id_buffer = indexToOriginalId;
        } else {
            temp_id_buffer.resize(num_vertices);
        }
        
        MPI_Bcast(temp_id_buffer.data(), num_vertices, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (rank != 0) {
            indexToOriginalId = temp_id_buffer;
        }
    }

    // SSSPGraph partitionMapitioning
    auto partitionMapition_start = high_resolution_clock::now();
    vector<idx_t> partitionMap(num_vertices);
    if (num_vertices > 0) {
        if (rank == 0) {
            cout << "Building CSR format for graph partitionMapitioning..." << endl;
        }
        
        vector<idx_t> xadjacencyList(num_vertices + 1);
        vector<idx_t> adjacencyListncy;
        
        if (rank == 0) {
            int edge_pos = 0;
            for (int i = 0; i < num_vertices; i++) {
                xadjacencyList[i] = edge_pos;
                for (GraphNode* curr = g.adjacencyList[i]; curr; curr = curr->next) {
                    adjacencyListncy.push_back(curr->dest);
                    edge_pos++;
                }
            }
            xadjacencyList[num_vertices] = edge_pos;
        }
        
        int edge_count = 0;
        if (rank == 0) {
            edge_count = adjacencyListncy.size();
            cout << "Number of edges: " << edge_count << endl;
        }
        MPI_Bcast(&edge_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (rank != 0) {
            adjacencyListncy.resize(edge_count);
        }
        
        MPI_Bcast(xadjacencyList.data(), num_vertices + 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (edge_count > 0) {
            MPI_Bcast(adjacencyListncy.data(), edge_count, MPI_INT, 0, MPI_COMM_WORLD);
        }
        
        if (rank != 0 && edge_count > 0) {
            for (int i = 0; i < num_vertices; i++) {
                for (int j = xadjacencyList[i]; j < xadjacencyList[i+1]; j++) {
                    g.addEdge(i, adjacencyListncy[j], 1);
                }
            }
        }
        
        if (rank == 0) {
            cout << "Partitioning graph with METIS..." << endl;
        }
        
        if (num_vertices > 0 && edge_count > 0) {
            idx_t ncon = 1;
            idx_t npartitionMaps = size;
            idx_t objval;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions(options);
            options[METIS_OPTION_CONTIG] = 0;
            
            int result = METIS_PartSSSPGraphKway(
                &num_vertices, &ncon, xadjacencyList.data(), adjacencyListncy.data(),
                NULL, NULL, NULL, &npartitionMaps,
                NULL, NULL, options, &objval, partitionMap.data()
            );
            
            if (result != METIS_OK) {
                cerr << "METIS partitionMapitioning failed on rank " << rank << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        } else {
            for (int i = 0; i < num_vertices; i++) {
                partitionMap[i] = i % size;
            }
        }
    }
    
    auto partitionMapition_end = high_resolution_clock::now();
    partitionMapition_time = duration_cast<milliseconds>(partitionMapition_end - partitionMapition_start).count();
    if (rank == 0) {
        cout << "Partitioning time: " << partitionMapition_time << " ms" << endl;
    }

    // Source node setup
    int src_id = -1, src = -1;
    if (rank == 0) {
        cout << "Enter source node ID: ";
        cin >> src_id;
        
        if (originalIdToIndex.find(src_id) != originalIdToIndex.end()) {
            src = originalIdToIndex[src_id];
            cout << "Source node " << src_id << " mapped to internal index " << src << endl;
        } else {
            cout << "Source node not found in graph." << endl;
        }
    }
    MPI_Bcast(&src, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (src == -1) {
        MPI_Finalize();
        return 0;
    }

    // Initialize distanceFromSourceances
    vector<int> distanceFromSource(num_vertices, INF);
    vector<int> predecessor(num_vertices, -1);
    if (src >= 0 && src < num_vertices && partitionMap[src] == rank) {
        distanceFromSource[src] = 0;
        predecessor[src] = -1;
    }

    // Initial SSSP computation
    if (num_vertices > 0) {
        if (rank == 0) {
            cout << "\n===== Initial SSSP Computation =====" << endl;
        }
        runDistributedDijkstra(g, distanceFromSource, predecessor, partitionMap, rank, size, total_sssp_time);
    }

    // Dynamic updates
    int num_updates = 0;
    if (rank == 0) {
        cout << "\nEnter number of updates: ";
        cin >> num_updates;
    }
    MPI_Bcast(&num_updates, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (num_vertices > 0) {
        srand(time(0) + rank * 17);
        
        for (int i = 0; i < num_updates; i++) {
            int update_type = 0, u = -1, v = -1;
            
            if (rank == 0) {
                cout << "\n===== Update " << i+1 << " =====" << endl;
            }
            
            auto update_start = high_resolution_clock::now();
            
            if (rank == 0) {
                update_type = rand() % 2;
                u = rand() % num_vertices;
                do {
                    v = rand() % num_vertices;
                } while (v == u);
                
                cout << (update_type ? "Delete" : "Insert") 
                     << " edge " << indexToOriginalId[u] << "->" << indexToOriginalId[v] << endl;
            }
            
            MPI_Bcast(&update_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&u, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&v, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            if (update_type == 0) {
                g.addEdge(u, v, 1);
            } else {
                g.removeEdge(u, v);
            }
            
            auto update_end = high_resolution_clock::now();
            update_time += duration_cast<milliseconds>(update_end - update_start).count();
            
            runDistributedDijkstra(g, distanceFromSource, predecessor, partitionMap, rank, size, total_sssp_time);
        }
    }

    writeProcessOutputToFile(distanceFromSource, predecessor, indexToOriginalId, src, rank);
    
    auto overall_end_time = high_resolution_clock::now();
    auto overall_duration = duration_cast<milliseconds>(overall_end_time - overall_start_time);
    
    if (rank == 0) {
        cout << "\n===== Performance Summary =====" << endl;
        cout << "Total execution time: " << overall_duration.count() << " milliseconds" << endl;
        cout << "  SSSPGraph reading time: " << graph_read_time << " milliseconds" << endl;
        cout << "  Partitioning time: " << partitionMapition_time << " milliseconds" << endl;
        cout << "  Dijkstra computation time: " << total_sssp_time << " milliseconds" << endl;
        cout << "  SSSPGraph update time: " << update_time << " milliseconds" << endl;
        cout << "  Other overhead: " << overall_duration.count() - graph_read_time - partitionMapition_time 
             - total_sssp_time - update_time << " milliseconds" << endl;
        cout << "Number of processes: " << size << endl;
        cout << "Number of vertices: " << num_vertices << endl;
        cout << "Number of dynamic updates: " << num_updates << endl;
    }

    g.clear();
    MPI_Finalize();
    return 0;
}

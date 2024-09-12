#include "Graph.h"
#include "MinHeap.cpp"
#include "MinHeap.h"

#include <iomanip>
#include <iostream>
#define INT_MAX 1000

template <typename T> struct greater {
  bool operator()(const T &x, const T &y) const { return x > y; }
};

template <typename T>
typename std::vector<T>::iterator find(std::vector<T> &vec, const T &value) {
  for (auto it = vec.begin(); it != vec.end(); ++it) {
    if (*it == value) {
      return it;
    }
  }
  return vec.end();
}

template <typename T> void Graph<T>::insert_vertex(const Vertex<T> &ver) {
  if (get_vertex_index(ver) == -1) {
    vertices.push_back(ver); // insert the vertex to the array of vertices
    std::vector<Edge> tmp;
    edges.push_back(tmp); // insert empty vector to the edges
  }
}

template <typename T>
int Graph<T>::get_vertex_index(const Vertex<T> &ver) const {
  for (int i = 0; i < vertices.size(); i++) {
    if (vertices[i].getOAirport() == ver.getOAirport()) {
      return i;
    }
  }

  return -1;
}

template <typename T>
void Graph<T>::add_edge(const Vertex<T> &ver1, const Vertex<T> &ver2,
                        int distance, int cost) {
  int i1 = get_vertex_index(ver1);
  int i2 = get_vertex_index(ver2);

  // Check if vertices exist, if not, add them
  if (i1 == -1) {
    insert_vertex(ver1);
    i1 = get_vertex_index(ver1);
  }
  if (i2 == -1) {
    insert_vertex(ver2);
    i2 = get_vertex_index(ver2);
  }

  if (i1 == -1 || i2 == -1) {
    throw std::string("Add_edge: incorrect vertices");
  }

  Edge v(i1, i2, distance, cost);

  edges[i1].push_back(v);
  if (i1 != i2) {
    edges[i2].push_back(v); // v2 should be v
  }
}

template <typename T> void Graph<T>::print() const {
  for (int i = 0; i < vertices.size(); i++) {
    std::cout << "{ " << vertices[i].getOAirport() << ": ";
    for (int j = 0; j < edges[i].size(); j++) {
      std::cout << '{' << vertices[edges[i][j].dest].getOAirport() << ", ";
      std::cout << edges[i][j].distance << "} ";
    }
    std::cout << " }\n";
  }
}

template <typename T> void reverse(std::vector<T> &vec) {
  int n = vec.size();
  for (int i = 0; i < n / 2; ++i) {
    std::swap(vec[i], vec[n - i - 1]);
  }
}

template <typename T> void Graph<T>::return_clean_visited() { clean_visited(); }

template <typename T> void Graph<T>::DFS(Vertex<T> &ver) {
  clean_visited();
  DFS_helper(ver);
  clean_visited();
}

template <typename T> void print_queue(std::vector<Vertex<T>> q) {
  while (!q.empty()) {
    std::cout << q.front().getOAirport() << " ";
    q.erase(q.begin());
  }
  std::cout << std::endl;
}

template <typename T> void Graph<T>::BFS(Vertex<T> &ver) {
  clean_visited();

  int i = get_vertex_index(ver);
  if (i == -1) {
    throw std::string("BFS: Vertex is not in the graph");
  }
  std::vector<int> q;
  q.push_back(i);
  vertices[i].setVisited(true);

  while (!q.empty()) {
    int i = q.front();
    std::cout << vertices[i].getOAirport() << ' ';
    for (int j = 0; j < edges[i].size(); j++) {
      int adjacent_ver = edges[i][j].dest;
      if (vertices[adjacent_ver].getVisited() == false) {
        vertices[adjacent_ver].setVisited(true);
        q.push_back(adjacent_ver);
      }
    }
    q.erase(q.begin());
  }

  clean_visited();
}

template <typename T> void Graph<T>::clean_visited() {
  for (Vertex<T> &v : vertices) {
    v.setVisited(false);
  }
}

template <typename T> void Graph<T>::DFS_helper(Vertex<T> &ver) {
  int i = get_vertex_index(ver);
  if (i == -1) {
    throw std::string("DFS: Vertex is not in the graph");
  }
  vertices[i].setVisited(true); // visit the current vertex
  std::cout << vertices[i].getOAirport()
            << ' '; // print the data of the current vertex
  for (int j = 0; j < edges[i].size(); j++) {
    Vertex<int> adjacent_ver = vertices[edges[i][j].dest];
    if (adjacent_ver.getVisited() == false) {
      DFS_helper(adjacent_ver);
    }
  }
}

//task 2
template <typename T>
void Graph<T>::dijkstra_shortest_path(const Vertex<T> &src,const Vertex<T> &dest) {
  int i_src = get_vertex_index(src);
  int i_dest = get_vertex_index(dest);

  if (i_src == -1 || i_dest == -1) {
    throw std::string("Shortest path: incorrect vertices");
  }

  clean_visited();
  std::vector<int> distances(vertices.size()); // represents the distances from
  // source to all other vertices
  std::vector<int> previous(
      vertices.size(), -1); // represents the previous vertex for each vertex

  // set initial distances
  for (int i = 0; i < distances.size(); i++) {
    distances[i] = (i == i_src) ? 0 : INT_MAX;
  }

  MinHeap<Edge> heap;
  int vertices_visited = 0;
  int cur_ver = i_src;

  while (vertices_visited < vertices.size()) {
    int i = cur_ver;
    // check the neighbours of the current node
    for (int j = 0; j < edges[i].size(); j++) {
      int i_adjacent_ver = edges[i][j].dest;
      if (vertices[i_adjacent_ver].getVisited() == false) {
        heap.insert(edges[i][j]);
        int dist_from_source = distances[i] + edges[i][j].distance;
        if (dist_from_source < distances[i_adjacent_ver]) {
          distances[i_adjacent_ver] = dist_from_source;
          previous[i_adjacent_ver] = i;
        }
      }
    }
    Edge e = heap.delete_min();
    cur_ver = e.dest;
    vertices[i].setVisited(true);
    vertices_visited++;
  }

  clean_visited();

  // Construct the shortest path
  std::vector<int> path;
  int cur_vertex = i_dest;
  while (cur_vertex != -1) {
    path.push_back(cur_vertex);
    cur_vertex = previous[cur_vertex];
  }

  for (size_t i = 0; i < path.size() / 2; ++i) {
    std::swap(path[i], path[path.size() - i - 1]);
  }

  // Print the shortest path
  if (distances[i_dest] == INT_MAX) {
    std::cout << "Shortest route from " << src.getOAirport() << " to "
              << dest.getOAirport() << ": None" << std::endl;
  } else {
    std::cout << "Shortest route from " << src.getOAirport() << " to "
              << dest.getOAirport() << ": ";
    for (size_t i = 0; i < path.size(); ++i) {
      std::cout << vertices[path[i]].getOAirport();
      if (i < path.size() - 1) {
        std::cout << " -> ";
      }
    }
    std::cout << ". The length is " << distances[i_dest] << ".";

    // Print the total cost
    int cost = 0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
      int index1 = path[i];
      int index2 = path[i + 1];
      for (const auto &edge : edges[index1]) {
        if (edge.dest == index2) {
          cost += edge.cost;
          break;
        }
      }
    }
    std::cout << " The total cost is " << cost << "." << std::endl;
  }
}

template <typename T>
std::vector<Edge> Graph<T>::get_edges(const Vertex<T> &ver) const {
  std::vector<Edge> incident_edges;

  // Find the index of the vertex
  int vertex_index = get_vertex_index(ver);
  if (vertex_index == -1) {
    std::cout << "get_edges: Vertex not found in the graph" << std::endl;
  } else {
    // Iterate over the edges incident to the vertex
    for (const Edge &edge : edges[vertex_index]) {
      incident_edges.push_back(edge);
    }
  }
  return incident_edges;
}

template <typename T> std::vector<Vertex<T>> Graph<T>::getVertices() const {
  return vertices;
}

template <typename T>
int Graph<T>::getEdgesCount(const Vertex<T> &vertex) const {
  int vertex_index = get_vertex_index(vertex);
  if (vertex_index == -1) {
    return 0;
  } else {
    return edges[vertex_index].size();
  }
}

template <typename T>
bool Graph<T>::edgeExists(const Vertex<T> &src, const Vertex<T> &dest) const {
  int src_index = get_vertex_index(src);
  int dest_index = get_vertex_index(dest);
  if (src_index == -1 || dest_index == -1) {
    return false;
  }
  for (const Edge &edge : edges[src_index]) {
    if (edge.dest == dest_index) {
      return true;
    }
  }
  return false;
}

template <typename T>
int Graph<T>::getEdgeCost(const Vertex<T> &src, const Vertex<T> dest) const {
  int src_index = get_vertex_index(src);
  int dest_index = get_vertex_index(dest);
  if (src_index == -1) {
    return -1;
  }
  for (const Edge &edge : edges[src_index]) {
    if (edge.dest == dest_index) {
      return edge.cost;
    }
  }
  return -1;
}

//task 6
template <typename T>
Graph<T> createUndirectedGraph(const Graph<T> &directedGraph) {
  Graph<T> undirectedGraph;

  // Iterate through all vertices in the directed graph
  for (const auto &vertex : directedGraph.getVertices()) {
    undirectedGraph.insert_vertex(vertex);
  }
  // Iterate through all edges in the directed graph
  for (const auto &vertex : directedGraph.getVertices()) {
    for (const auto &edge : directedGraph.get_edges(vertex)) {
      int j = edge.dest;
      // Check if there is only one directed edge between vertices i and j
      if (!undirectedGraph.edgeExists(vertex, directedGraph.getVertices()[j])) {
        undirectedGraph.add_edge(vertex, directedGraph.getVertices()[j],
                                 edge.distance, edge.cost);
      }
      if (!undirectedGraph.edgeExists(directedGraph.getVertices()[j], vertex)) {
        undirectedGraph.add_edge(directedGraph.getVertices()[j], vertex,
                                 edge.distance, edge.cost);
      }
    }
  }
  return undirectedGraph;
}

template <typename T> void print_path(const std::vector<int> &path) {
  if (path.empty()) {
    std::cout << "None" << std::endl;
  }
  for (size_t i = 0; i < path.size() - 1; ++i) {
    std::cout << path[i] << "-> ";
  }
  std::cout << path.back() << std::endl;
}

template <typename T> bool Graph<T>::isUndirected() const {
  for (int i = 0; i < edges.size(); ++i) {
    for (int j = 0; j < edges[i].size(); ++j) {
      int src = i;
      int dest = edges[i][j].dest;

      bool found = false;
      for (int k = 0; k < edges[dest].size(); ++k) {
        if (edges[dest][k].dest == src) {
          found = true;
          break;
        }
      }
      if (!found) {
        return false;
      }
    }
  }

  return true;
}

template <typename T>
std::vector<Edge> Graph<T>::return_edges(const Vertex<T> &ver) const {
  std::vector<Edge> edges_vector;
  int vertex_index = get_vertex_index(ver);
  if (vertex_index == -1) {
    std::cout << "Not found\n";
  } else {
    for (const Edge &edge : edges[vertex_index]) {
      edges_vector.push_back(edge);
    }
  }
  return edges_vector;
}

template <typename T>
const std::vector<Vertex<T>> &Graph<T>::return_vertices() const {
  return vertices;
}

template <typename T> struct ShortestPathResult {
  std::vector<int> path;
  int length;
};

// task 3
template <typename T>
std::vector<int> Graph<T>::findVertWithST(const std::string &state) {
  std::vector<int> stateVerts;

  for (size_t i = 0; i < vertices.size(); ++i) {
    int pos = vertices[i].getOCity().find_last_of(' ');
    std::string temp = vertices[i].getOCity().substr(pos + 1, 2);
    if (temp == state) {
      stateVerts.push_back(i);
    }
  }
  return stateVerts;
}

template <typename T>
void Graph<T>::shortest_paths_to_state(const Vertex<T> &src,
                                       const std::string &state) {
  int i_src = get_vertex_index(src);

  std::vector<int> stateVerts = findVertWithST(state);

  if (i_src == -1) {
    throw std::string("Incorrect source");
  }

  clean_visited();
  std::vector<int> distances(vertices.size());
  std::vector<int> previous(vertices.size(), -1);

  for (int i = 0; i < distances.size(); i++) {
    distances[i] = (i == i_src) ? 0 : INT_MAX;
  }

  MinHeap<Edge> heap;
  int vertices_visited = 0;
  int cur_ver = i_src;

  while (vertices_visited < vertices.size()) {
    int i = cur_ver;
    // check the neighbours of the current node
    for (int j = 0; j < edges[i].size(); j++) {
      int i_adjacent_ver = edges[i][j].dest;
      if (vertices[i_adjacent_ver].getVisited() == false) {
        heap.insert(edges[i][j]);
        int dist_from_source = distances[i] + edges[i][j].distance;
        if (dist_from_source < distances[i_adjacent_ver]) {
          distances[i_adjacent_ver] = dist_from_source;
          previous[i_adjacent_ver] = i;
        }
      }
    }
    Edge e = heap.delete_min();
    cur_ver = e.dest;
    vertices[i].setVisited(true);
    vertices_visited++;
  }

  clean_visited();

  std::cout << "Shortest route from " << src.getOAirport() << " to " << state
            << " State are: " << std::endl;
  std::cout << "Path                Length     Cost" << std::endl;
  // Construct the shortest path
  for (int i = 0; i < stateVerts.size(); i++) {
    std::vector<int> path;
    int cur_vertex = stateVerts[i]; // start from the destination vertex
    while (cur_vertex != -1) {
      path.push_back(cur_vertex);
      cur_vertex = previous[cur_vertex];
    }
    // path.push_back(i_src); // add the source vertex to the path

    // Print the shortest path
    if (distances[stateVerts[i]] == INT_MAX) {
      // nothing
    } else {
      for (size_t k = path.size() - 1; k > 0; --k) {
        std::cout << vertices[path[k]].getOAirport();
        if (k >= 1) {
          std::cout << " -> ";
        }
      }
      std::cout << vertices[path[0]].getOAirport();
      std::cout << "           " << distances[stateVerts[i]] << "       ";

      // Print the total cost
      int cost = 0;
      for (size_t p = path.size() - 1; p > 0; --p) {
        int index1 = path[p];
        int index2 = path[p - 1];
        for (const auto &edge : edges[index1]) {
          if (edge.dest == index2) {
            cost += -edge.cost;
            break;
          }
        }
      }
      std::cout << -cost << std::endl;
    }
  }
}

// task 4
template <typename T>
void DFS_shortest_path(const Graph<T> &graph, const Vertex<T> &src,
                       const Vertex<T> &dest, int num_stops,
                       std::vector<int> &path, std::vector<int> &shortest_path,
                       int &shortest_length, int &shortest_cost, int stops,
                       int &current_length, int &current_cost) {
  int src_index = graph.get_vertex_index(src);
  int dest_index = graph.get_vertex_index(dest);

  // Check if we reached the destination
  if (src_index == dest_index && stops <= num_stops) {
    // If we have a better solution, update the shortest path and its length and
    // cost
    if (path.size() < shortest_length) {
      shortest_path = path;
      shortest_length = current_length;
      shortest_cost = current_cost;
    }
    return;
  }

  // Iterate over all possible flights from the current source vertex
  for (const Edge &edge : graph.return_edges(src)) {
    // Check if the destination vertex is already in the path
    if (find(path, edge.dest) == path.end()) {
      // Add the current vertex to the path
      path.push_back(edge.dest);
      // Update the current length and cost
      current_length += edge.distance;
      current_cost += edge.cost;
      // Recursive call to explore further
      DFS_shortest_path(graph, graph.return_vertices()[edge.dest], dest,
                        num_stops, path, shortest_path, shortest_length,
                        shortest_cost, stops + 1, current_length, current_cost);
      // Backtrack: remove the last vertex from the path
      path.pop_back();
      current_length -= edge.distance;
      current_cost -= edge.cost;
    }
  }
}

template <typename T>
void shortest_path_with_stops(Graph<T> &graph, const Vertex<T> &src,
                              const Vertex<T> &dest, int num_stops) {
  std::vector<int> path;
  std::vector<int> shortest_path;
  int shortest_length = INT_MAX;
  int shortest_cost = INT_MAX;

  // Add the source vertex to the path
  path.push_back(graph.get_vertex_index(src));
  int current_length = 0;
  int current_cost = 0;

  // Call the recursive function to find the shortest path with the specified
  // number of stops
  DFS_shortest_path(graph, src, dest, num_stops, path, shortest_path,
                    shortest_length, shortest_cost, 0, current_length,
                    current_cost);

  // Print the result
  if (shortest_length == INT_MAX) {
    std::cout << "Shortest route from " << src.getOAirport() << " to "
              << dest.getOAirport() << " with " << num_stops << " stops: None"
              << std::endl;
  } else {
    std::cout << "Shortest route from " << src.getOAirport() << " to "
              << dest.getOAirport() << " with " << num_stops << " stops: ";
    for (size_t i = 0; i < shortest_path.size(); ++i) {
      std::cout << graph.getVertices().at(shortest_path[i]).getOAirport();
      if (i < shortest_path.size() - 1) {
        std::cout << "->";
      }
    }
    std::cout << " The length is " << shortest_length - 1 << "."
              << " The cost is " << shortest_cost << "." << std::endl;
  }
}

// task 5
template<typename T>
void print_num_edges(const Graph<T>& graph, Vertex<T> vertex) {
  std::vector<Edge> edge = graph.get_edges(vertex);

  std::cout<< vertex.getOAirport()<< "              "<< edge.size()<<std::endl;
}

//task 7
template <typename T>
void Graph<T>::primMST(std::vector<Edge> &mst, int &totalCost) const {
  int numVertices = vertices.size();
  std::vector<bool> visitedV(numVertices, false);
  mst.clear();
  totalCost = 0;
  if (numVertices == 0) {
    std::cout << "Graph is empty\n";
    return;
  }
  visitedV[0] = true;
  while (mst.size() < numVertices - 1) {
    int minCost = INT_MAX;
    Edge minEdge;
    bool found = false;
    for (int i = 0; i < numVertices; ++i) {
      if (visitedV[i]) {
        for (const Edge &edge : edges[i]) {
          if (!visitedV[edge.dest] && edge.cost < minCost) {
            minCost = edge.cost;
            minEdge = edge;
            found = true;
          }
        }
      }
    }
    if (!found) {
      return;
    }
    mst.push_back(minEdge);
    totalCost += minEdge.cost;
    visitedV[minEdge.dest] = true;
  }
}

template <typename T>
void bubbleSort(std::vector<T> &vect, bool (*comp)(const T &, const T &)) {
  int n = vect.size();
  for (int i = 0; i < n - 1; ++i) {
    for (int j = 0; j < n - i - 1; ++j) {
      if (comp(vect[j], vect[j + 1])) {
        std::swap(vect[j], vect[j + 1]);
      }
    }
  }
}
// comparison for kruskalMST
template <typename T>
bool Graph<T>::compareEdges(const Edge &a, const Edge &b) {
  return a.cost < b.cost;
}

//task 8
template <typename T>
void Graph<T>::kruskalMST(std::vector<Edge> &mst, int &totalCost) const {
  int verticesSize = vertices.size();
  mst.clear();
  totalCost = 0;
  // create vector to store all edges and sort
  std::vector<Edge> sortedEdges;
  for (const auto &adjList : edges) {
    for (const auto &edge : adjList) {
      sortedEdges.push_back(edge);
    }
  }
  bubbleSort(sortedEdges, compareEdges);

  std::vector<int> parent(verticesSize);
  for (int i = 0; i < verticesSize; ++i) {
    parent[i] = i;
  }

  auto find = [&](int k) {
    while (k != parent[k]) {
      parent[k] = parent[parent[k]];
      k = parent[k];
    }
    return k;
  };
  for (const auto &edge : sortedEdges) {
    int sParent = find(edge.src);
    int dParent = find(edge.dest);
    if (sParent != dParent) {
      mst.push_back(edge);
      totalCost += edge.cost;
      parent[sParent] = dParent;
    }
  }
}

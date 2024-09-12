#ifndef GRAPH_H
#define GRAPH_H
#include "Vertex.h"
#include "Edge.h"
#include <vector>
#include <string>


template <typename T>
class Graph {
public:
    Graph() {}



    

    void insert_vertex(const Vertex<T>& ver);
    void add_edge(const Vertex<T>& ver1, const Vertex<T>& ver2, int distance, int cost); //connect ver1 with ver2
    void return_clean_visited();
    void print() const;
    void countDirectFlights();
    void DFS(Vertex<T>& ver);
    void BFS(Vertex<T>& ver);
    void dijkstra_shortest_path(const Vertex<T>& src, const Vertex<T>& dest); 

    int get_vertex_index(const Vertex<T>& ver) const;
    std::vector<Edge> return_edges(const Vertex<T>& ver) const;
    bool isUndirected() const;

    const std::vector<Vertex<T>>& return_vertices() const;
    std::vector<Vertex<T>> getVertices() const;
    std::vector<Edge> get_edges(const Vertex<T>& ver) const;
    int getEdgesCount(const Vertex<T>& vertex)const;
    bool edgeExists(const Vertex<T>& src, const Vertex<T>& dest) const;
    int getEdgeCost(const Vertex<T>& src, const Vertex<T> dest) const;
    void kruskalMST(std::vector<Edge>& mst, int& totalCost ) const;
    void primMST(std::vector<Edge>& mst, int& totalCost ) const;
    void shortest_paths_to_state(const Vertex<T>& src, const std::string& state);
    static bool compareEdges(const Edge& a, const Edge& b);

std::string getVertexAirport(int index) const{
  if(index >=0 && index < vertices.size()){
    return vertices[index].getOAirport();
  }
  return "";
}

private:
  std::vector<Vertex<T>> vertices; //nodes
  std::vector<std::vector<Edge>> edges; //connections
  std::vector<int> directFlightConnections;
  void clean_visited();

  void DFS_helper(Vertex<T>& ver);
  std::vector<int> findVertWithST(const std::string& state);

};

#endif

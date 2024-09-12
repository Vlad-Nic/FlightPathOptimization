#include "Graph.h"
#include "Vertex.h"
#include "Graph.cpp"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;

int main() {
  ifstream inFile("airports.csv");

  string orAirport, dAirport, orCity, dCity;
  int distance, cost;

  string title1;
  string row;

  getline(inFile, title1);
  
  Graph<string> airportGraph;

    //create graph based on the file
  while (getline(inFile, row)) {
    stringstream ss(row);
    string field;
    string filler1, filler2;

    getline(ss, orAirport, ',');
    getline(ss, dAirport, ',');
    
    getline(ss, filler1, '"');
    getline(ss, orCity, '"');
    getline(ss, filler1, ',');

    getline(ss, filler2, '"');
    getline(ss, dCity, '"'); 
    getline(ss, filler2, ',');
    
    getline(ss, field, ',');
    distance = atoi(field.c_str());
    
    getline(ss, field, ',');
    cost = atoi(field.c_str());

    Vertex<string> originVertex(orAirport, orCity);
    Vertex<string> destinationVertex(dAirport, dCity);

    // Add vertices to the graph only if they don't already exist
    if (airportGraph.get_vertex_index(originVertex) == -1) {
      airportGraph.insert_vertex(originVertex);
    }

    if (airportGraph.get_vertex_index(destinationVertex) == -1) {
      airportGraph.insert_vertex(destinationVertex);
    }

    // Add an edge between the origin and destination vertices
    airportGraph.add_edge(originVertex, destinationVertex, distance, cost);
  }

  inFile.close();

  // airportGraph.print();

  
  // test task 2
  Vertex<string> source("ROC","Rochester, NY");
  Vertex<string> dest("MIA", "Miami, FL");
  airportGraph.dijkstra_shortest_path(source, dest);

  std::cout<<std::endl;

  Vertex<string> source_("AUS","Austin, TX");
  Vertex<string> dest_("ORD", "Chicago, IL");
  airportGraph.dijkstra_shortest_path(source_, dest_);

  std::cout<<std::endl;
  
  //test task 3 
  Vertex<std::string> src("ATL","Atlanta");
  string st = "FL";
  airportGraph.shortest_paths_to_state(src,st);

  std::cout<<std::endl;
  
  //test 4
  
  Vertex<string> source2("BOS","BOSTON, MA");
  Vertex<string> dest2("DTW", "Detroit, MI");
  shortest_path_with_stops(airportGraph, source2, dest2, 1);
  std::cout<<std::endl;
  
  //test 5
  std::cout<< "Airport        Connections"<<std::endl;
  print_num_edges(airportGraph, source2);
  print_num_edges(airportGraph, source);
  print_num_edges(airportGraph, dest);

  std::cout<<std::endl;
  
  //test task 6
  Graph<string> undirected = createUndirectedGraph(airportGraph);
  // undirected.print();
  //test task 7
  std::vector<Edge> mstEdges;
  int totalCost;
  undirected.primMST(mstEdges, totalCost);
  std::cout << "Minimal Spanning Tree:\n";
  for(const auto& edge : mstEdges){
     std:: cout << "Edge: " << undirected.getVertexAirport(edge.src) << " - " << undirected.getVertexAirport(edge.dest) << " " << "Weight: " << edge.cost << std::endl;
  }
  std::cout << "Total Cost of MST: " << totalCost << std::endl<<std::endl;

  //test 8 
  std::vector<Edge> kruskalEdges;
  int kruskalCost;
  undirected.kruskalMST(kruskalEdges, kruskalCost);
  std::cout << "Minimal Spanning Tree:\n";
  for(const auto& edge : kruskalEdges){
    std:: cout << "Edge: " << undirected.getVertexAirport(edge.src) << " - " << undirected.getVertexAirport(edge.dest) << " " << "Weight: " << edge.cost << std::endl;
  }
  std::cout << "Total Cost of MST: " << kruskalCost << std::endl;
  
  return 0;
}







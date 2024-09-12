# FlightPathOptimization
This project aims to develop a comprehensive system for modeling and optimizing airport connectivity and flight routes. 
Tasks
1.	Airport Connectivity Map: The system constructs a graph where airports are nodes, and direct flights are edges. This map serves as the foundation for analyzing the connectivity of the airport network, enabling tasks such as finding the shortest path between airports.
2.	Flight Route Optimization: Building on the connectivity map, the system employs graph algorithms like DFS, BFS, Dijkstra's, Prim’s and Kruskal’s to find the most efficient routes between airports. This feature helps in identifying optimal routes for passengers and cargo, potentially reducing travel time and costs.
3.	Efficiency Check/ Evaluation of the Algorithms:

Functionality
1.	Read the information from the dataset (a csv file) and create a weighted directed graph G. Note that you need to consider two weights for each edge. One is the Distance, and the other is Cost.
2.	Find the shortest path between the given origin airport and destination airport. The algorithm should output both the path and the length of the path. The algorithm should provide the appropriate message if such path doesn’t exist.
3.	Find all shortest paths between the given origin airport and all the airports in the destination state. The algorithm should output all the paths and their lengths. The algorithm should provide the appropriate message if no such paths exist.
4.	Find the shortest path between the given origin airport and destination airport with a given number of stops. The algorithm should provide the appropriate message if such a path doesn’t exist.
5.	Count and display total direct flight connections to each airport. You should consider both outbound and inbound flights. For instance, if you can directly fly to Tampa airport only from Miami, Orlando, and Atlanta, the inbound count for Tampa airport would be 3. If you can directly fly from Tampa airport only to New York, the outbound count for Tampa airport is 1. The total number of direct flight connections for Tampa airport is 4. The list of airports should be sorted based on the total direct flight connections count, starting with the airports having the highest number of direct flight connections.
6.	Create an undirected graph G_u from the original directed graph G using the following rules:
a.	For each pair of vertices u and v, if there is only one directed edge(either (u,v) or (v,u)) between them, you keep that single edge with its corresponding cost as an undirected weighted edge. You can ignore the distance on that edge.
b.	For each pair of vertices u and v, if there are two directed edges (u,v) and (v, u) between them, you keep the one with the minimum cost value as an undirected weighted edge. You can ignore the distance on that edge.
7.	Generate a Minimal Spanning Tree utilizing Prim’s algorithm on G_u that you created in the previous step. The algorithm will output both the content of the constructed MST and its total cost. In this step, for each edge you need to consider the cost as weight to minimize the total cost. In the event of a disconnected graph, the algorithm will appropriately notify that an MST cannot be formed. Note: A connected graph is defined as one where there exists a path between every pair of
vertices.
8.	Generate a Minimal Spanning Tree using Kruskal’s algorithm on G_u that you created in the previous step. The algorithm will output both the content of the constructed MST and its total cost. In this step, for each edge you need to consider the cost as weight to minimize the total cost. If the graph is disconnected the algorithm should provide minimum spanning forest consisting of a minimum spanning tree for each connected component.

Classes Used

•	Vertex: Represents a vertex in the graph. It contains information about the airport code and city.

•	Edge: Represents an edge between two vertices in the graph. It contains information about the source vertex, destination vertex, distance, and cost.

•	Graph: Implements the graph data structure using an adjacency list representation. It includes methods for inserting vertices, adding edges, performing traversals (DFS, BFS), finding shortest paths, and more.

•	MinHeap: Implements a min-heap data structure, used primarily in Dijkstra's algorithm for maintaining the priority queue of vertices based on their distances.

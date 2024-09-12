#ifndef VERTEX_H
#define VERTEX_H

template <typename T>
class Vertex {
public:
    Vertex(const T& oAirport = T(), const T& oCity = T()) : oAirport(oAirport), oCity(oCity), visited(false) {};
    const T& getOAirport() const {return oAirport; }
    const T& getOCity() const {return oCity; }
    bool getVisited() const {return visited; }
    void setVisited(bool v) { visited = v; }

private:
  T oAirport;
  T oCity;
  bool visited;
};

#endif
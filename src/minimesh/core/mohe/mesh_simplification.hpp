#pragma once

//
// Mesh_modifier_for_simplification.hpp
//
// Some functionality for modifying meshes.
//
// Author: Shayan Hoshyari
//


#include <string>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{

// Forward declaration
class Simplifier;

// Structure to represent a quadric (4x4 matrix for error metric)
struct Quadric 
{
    Eigen::Matrix4d Q;
    
    Quadric();
    Quadric(const Eigen::Vector4d& plane);
    Quadric operator+(const Quadric& other) const;
    double computeError(const Eigen::Vector3d& v) const;
    Eigen::Vector3d findOptimalPosition(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const;
};

// Structure to represent an edge collapse operation
struct EdgeCollapse 
{
    int he_index;           // Half-edge to collapse
    int v1, v2;             // Vertices of the edge (v1 will be kept, v2 will be removed)
    double cost;            // Cost of the collapse
    Eigen::Vector3d new_pos; // New position for the remaining vertex
    
    // For priority queue (min-heap based on cost)
    bool operator>(const EdgeCollapse& other) const;
};


class Mesh_modifier_for_simplification
{
public:
	// Trivial constructor
	Mesh_modifier_for_simplification(Mesh_connectivity & mesh_in): _m(mesh_in) {}

	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	//
	// Given two vertices, this function return the index of the half-edge going from v0 to v1.
	// Returns mesh::invalid_index if no half-edge exists between the two vertices.
	//
	int get_halfedge_between_vertices(Mesh_connectivity& mesh, const int v0, const int v1);

	//
	// Flip an edge in a mesh
	// Input: The mesh, and the index of a half-edge on the edge we wish to flip
	// Return true if operation successful, and false if operation not possible
	//
	// Assumption: mesh is all triangles
	//
	// NOTE: To see how this method works, take a look at edge-flip.svg
	//
	bool flip_edge(const int he_index);


private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;
};


class Simplifier
{
public:
	// Trivial constructor
	Simplifier(Mesh_connectivity & mesh_in): _m(mesh_in), _proxy_mesh() {}

	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }
	
	// Get the proxy mesh
	Mesh_connectivity & proxy_mesh() { return _proxy_mesh; }
	const Mesh_connectivity & proxy_mesh() const { return _proxy_mesh; }

	bool simplify_try_before_commit(int target_faces); 
  bool simplify_test_ahead(int target_faces);
  bool simplify_once();

private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;
	// proxy mesh for testing operations
	Mesh_connectivity _proxy_mesh;
	
	// Helper methods for QSlim algorithm
	void initializeVertexQuadrics(std::map<int, Quadric>& vertex_quadrics);
	
	void buildInitialEdgeQueue(const std::map<int, Quadric>& vertex_quadrics,
	                          std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, std::greater<EdgeCollapse>>& edge_queue,
	                          std::set<std::pair<int,int>>& processed_edges);
	
	EdgeCollapse computeEdgeCollapse(int v1, int v2, int he_idx, const std::map<int, Quadric>& vertex_quadrics);

	void checkVertexAssociatedHalfEdges(Mesh_connectivity& mesh, int halfedge_idx);

	bool isValidCollapse(const EdgeCollapse& collapse);
	
	bool performEdgeCollapse(Mesh_connectivity& mesh, EdgeCollapse& collapse, std::map<int, Quadric>& vertex_quadrics);
	// bool performEdgeCollapseOnProxy(EdgeCollapse& collapse, std::map<int, Quadric>& vertex_quadrics);
	
	void initializeProxyMesh();
	void copyMeshToProxy();
	void copyProxyToMesh();
	bool checkProxyManifold(Mesh_connectivity& mesh);

	bool collapseEdgeTopology(Mesh_connectivity & mesh, int he_idx, int v_keep, int v_remove);
	// bool collapseEdgeTopologyOnProxy(int he_idx, int v_keep, int v_remove);
	
	void updateEdgeQueue(int vertex_idx, const std::map<int, Quadric>& vertex_quadrics,
	                    std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, std::greater<EdgeCollapse>>& edge_queue,
	                    std::set<std::pair<int,int>>& processed_edges);
	
	bool wouldCreateNonManifold(int v1, int v2);
	int get_halfedge_between_vertices(Mesh_connectivity& mesh, const int v0, const int v1);

};


} // end of minimesh
}

# Mesh Connectivity (minimesh::mohecore)

The `Mesh_connectivity` class implements a **half-edge data structure** for polygonal meshes.  
It provides **iterators** for vertices, half-edges, faces, and rings, as well as utilities for modifying, validating, and copying meshes.

---

## 1. Iterators

Iterators are lightweight handles to mesh elements.

### 🔹 Vertex Iterator
| Function | Description |
|----------|-------------|
| `Vertex_data & data()` | Access vertex data (coordinates, connectivity). |
| `int index()` | Get the vertex index. |
| `bool is_active()` | Check if the vertex is alive. |
| `void deactivate()` | Mark vertex as deleted. |
| `Half_edge_iterator half_edge()` | Get one outgoing half-edge. |
| `const Eigen::Vector3d & xyz()` | Get vertex position. |
| `bool is_equal(Vertex_iterator)` | Compare with another iterator. |

---

### 🔹 Half-edge Iterator
| Function | Description |
|----------|-------------|
| `Half_edge_data & data()` | Access half-edge data. |
| `int index()` | Get half-edge index. |
| `bool is_active()` | Check if the half-edge is alive. |
| `void deactivate()` | Mark half-edge as deleted. |
| `Half_edge_iterator next()` | Get the next half-edge in the face loop. |
| `Half_edge_iterator prev()` | Get the previous half-edge. |
| `Half_edge_iterator twin()` | Get the twin half-edge (opposite direction). |
| `Vertex_iterator origin()` | Get the start vertex. |
| `Vertex_iterator dest()` | Get the end vertex. |
| `Face_iterator face()` | Get the incident face. |
| `bool is_equal(Half_edge_iterator)` | Compare with another iterator. |

---

### 🔹 Face Iterator
| Function | Description |
|----------|-------------|
| `Face_data & data()` | Access face data. |
| `int index()` | Get face index. |
| `Half_edge_iterator half_edge()` | Get one boundary half-edge. |
| `int n_vertices()` | Count vertices of the face. |
| `bool is_active()` | Check if face is alive. |
| `void deactivate()` | Mark face as deleted. |
| `bool is_equal(Face_iterator)` | Compare with another iterator. |

---

### 🔹 Vertex Ring Iterator
Iterates over all half-edges **incoming to a vertex** (the “star” around a vertex).

| Function | Description |
|----------|-------------|
| `Half_edge_iterator half_edge()` | Current half-edge. |
| `bool reset_boundary()` | Reset to boundary edge if vertex lies on boundary. |
| `bool advance()` | Move to next half-edge around the vertex, returns false if looped. |

---

## 2. Mesh Access & Creation

| Function | Description |
|----------|-------------|
| `Vertex_iterator vertex_at(int vertex_id)` | Get vertex iterator. |
| `Vertex_ring_iterator vertex_ring_at(int vertex_id)` | Get vertex ring iterator. |
| `Face_iterator face_at(int face_id)` | Get face iterator. |
| `Half_edge_iterator half_edge_at(int he_id)` | Get half-edge iterator. |
| `Face_iterator hole()` | Get the special “hole” face (for boundaries). |

---

## 3. Counting Entities

| Function | Description |
|----------|-------------|
| `int n_active_vertices()` | Number of active vertices. |
| `int n_active_half_edges()` | Number of active half-edges. |
| `int n_active_faces()` | Number of active faces. |
| `int n_total_vertices()` | Total allocated vertices (active + inactive). |
| `int n_total_half_edges()` | Total allocated half-edges. |
| `int n_total_faces()` | Total allocated faces. |

---

## 4. Modifying Mesh

| Function | Description |
|----------|-------------|
| `Vertex_iterator add_vertex(bool allow_recycling = true)` | Add a new vertex. |
| `Half_edge_iterator add_half_edge(bool allow_recycling = true)` | Add a new half-edge. |
| `Face_iterator add_face(bool allow_recycling = true)` | Add a new face. |

---

## 5. Defragmentation

| Function | Description |
|----------|-------------|
| `void compute_defragmention_maps(Defragmentation_maps & m)` | Compute maps for compacting indices. |
| `void defragment(Defragmentation_maps & m, Mesh_connectivity & new_mesh)` | Build a new defragmented mesh. |
| `void defragment_in_place(Defragmentation_maps & m)` | Defragment the mesh in place. |

---

## 6. Building Mesh Connectivity

| Function | Description |
|----------|-------------|
| `void build_from_triangles(const std::vector<double>& xyz, const std::vector<int>& triangles)` | Build mesh from triangles. |
| `void build_from_polygons(const std::vector<double>& xyz, const std::vector<int>& polygon_verts, const std::vector<int>& polygon_adj)` | Build mesh from polygons. |

---

## 7. Validation

| Function | Description |
|----------|-------------|
| `bool check_sanity_slowly(bool verbose = false)` | Check mesh consistency (slow but thorough). |

---

## 8. Copying & Clearing

| Function | Description |
|----------|-------------|
| `void copy(const Mesh_connectivity & other)` | Copy another mesh. |
| `void swap(Mesh_connectivity & other)` | Swap contents with another mesh. |
| `void clear()` | Clear all data. |

---

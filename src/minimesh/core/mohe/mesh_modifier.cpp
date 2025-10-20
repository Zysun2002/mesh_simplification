#include <minimesh/core/mohe/mesh_modifier.hpp>
#include <minimesh/core/util/assert.hpp>

#include <iostream>
#include <map>
#include <vector>

#define M_PI 3.14159265358979323846


namespace minimesh
{
namespace mohecore
{


//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int Mesh_modifier::get_halfedge_between_vertices(const int v0, const int v1)
{
	// Get a ring iterator for v0
	Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(v0);

	int answer = mesh().invalid_index;

	// Loop over all half-edges that end at v0.
	do
	{
		// Make sure that the half-edge does end at v0
		assert(ring_iter.half_edge().dest().index() == v0);

		// If the half-edge also starts and v1, then it's twin
		// goes from v0 to v1. This would be the half-edge that
		// we were looking for
		if(ring_iter.half_edge().origin().index() == v1)
		{
			answer = ring_iter.half_edge().twin().index();
		}
	} while(ring_iter.advance());

	if(answer != mesh().invalid_index)
	{
		assert(mesh().half_edge_at(answer).origin().index() == v0);
		assert(mesh().half_edge_at(answer).dest().index() == v1);
	}

	return answer;
}


bool Mesh_modifier::flip_edge(const int he_index)
{
	//
	// Take a reference to all involved entities
	//

	// HALF-EDGES
	Mesh_connectivity::Half_edge_iterator he0 = mesh().half_edge_at(he_index);
	Mesh_connectivity::Half_edge_iterator he1 = he0.twin();

	// meshes on the boundary are not flippable
	if(he0.face().is_equal(mesh().hole()) || he1.face().is_equal(mesh().hole()))
	{
		return false;
	}

	Mesh_connectivity::Half_edge_iterator he2 = he0.next();
	Mesh_connectivity::Half_edge_iterator he3 = he2.next();
	Mesh_connectivity::Half_edge_iterator he4 = he1.next();
	Mesh_connectivity::Half_edge_iterator he5 = he4.next();

	// VERTICES
	Mesh_connectivity::Vertex_iterator v0 = he1.origin();
	Mesh_connectivity::Vertex_iterator v1 = he0.origin();
	Mesh_connectivity::Vertex_iterator v2 = he3.origin();
	Mesh_connectivity::Vertex_iterator v3 = he5.origin();

	// FACES
	Mesh_connectivity::Face_iterator f0 = he0.face();
	Mesh_connectivity::Face_iterator f1 = he1.face();

	//
	// Now modify the connectivity
	//

	// HALF-EDGES
	he0.data().next = he3.index();
	he0.data().prev = he4.index();
	he0.data().origin = v3.index();
	//
	he1.data().next = he5.index();
	he1.data().prev = he2.index();
	he1.data().origin = v2.index();
	//
	he2.data().next = he1.index();
	he2.data().prev = he5.index();
	he2.data().face = f1.index();
	//
	he3.data().next = he4.index();
	he3.data().prev = he0.index();
	//
	he4.data().next = he0.index();
	he4.data().prev = he3.index();
	he4.data().face = f0.index();
	//
	he5.data().next = he2.index();
	he5.data().prev = he1.index();

	// VERTICES
	v0.data().half_edge = he2.index();
	v1.data().half_edge = he4.index();
	v2.data().half_edge = he1.index();
	v3.data().half_edge = he0.index();

	// FACES
	f0.data().half_edge = he0.index();
	f1.data().half_edge = he1.index();

	// operation successful
	return true;
} // All done


bool Butterfly_subdivider::subdivide_loop()
{
    using namespace mohecore;

    std::map<int, Mesh_connectivity::Vertex_iterator> edge_midpoints;

    // create new vertices at midpoints
	for(int heid = 0; heid < _m.n_total_half_edges(); ++heid)
	{
		Mesh_connectivity::Half_edge_iterator he = _m.half_edge_at(heid);
		
		if(edge_midpoints.count(he.index()) > 0) continue;

		Mesh_connectivity::Half_edge_iterator twin = he.twin();

		Eigen::Vector3d v0 = he.origin().xyz();
		Eigen::Vector3d v1 = he.dest().xyz();

		Mesh_connectivity::Vertex_iterator v2 = he.next().dest();     
		Mesh_connectivity::Vertex_iterator v3 = twin.next().dest();   

		Mesh_connectivity::Vertex_iterator v4 = he.prev().twin().next().dest();    
		Mesh_connectivity::Vertex_iterator v5 = he.next().twin().next().dest(); 
		Mesh_connectivity::Vertex_iterator v6 = twin.prev().twin().next().dest();    
		Mesh_connectivity::Vertex_iterator v7 = twin.next().twin().next().dest(); 

		Eigen::Vector3d new_pos = 0.5*(v0 + v1)
								+ 0.125*(v2.xyz() + v3.xyz())
								- 0.0625*(v4.xyz() + v5.xyz() + v6.xyz() + v7.xyz());

		// register new vertex
		Mesh_connectivity::Vertex_iterator new_vertex = _m.add_vertex();
		new_vertex.data().xyz = new_pos;

		// update half edge
		edge_midpoints[he.index()] = new_vertex;
		edge_midpoints[twin.index()] = new_vertex;
	}

    // split faces to get new sub-faces
    std::vector<int> new_triangles;
    for(int fid = 0; fid < _m.n_total_faces(); ++fid)
    {
        Mesh_connectivity::Face_iterator face = _m.face_at(fid);
        if(!face.is_active()) continue;//

        Mesh_connectivity::Half_edge_iterator he0 = face.half_edge();
        Mesh_connectivity::Half_edge_iterator he1 = he0.next();
        Mesh_connectivity::Half_edge_iterator he2 = he0.prev();

        int v0 = he0.origin().index();
        int v1 = he1.origin().index();
        int v2 = he2.origin().index();

        int m0 = edge_midpoints[he0.index()].index();
        int m1 = edge_midpoints[he1.index()].index();
        int m2 = edge_midpoints[he2.index()].index();

        new_triangles.insert(new_triangles.end(), {v0,m0,m2, v1,m1,m0, v2,m2,m1, m0,m1,m2});
        face.deactivate();
    }

    // reconstruct mesh
    std::vector<double> xyz;
    xyz.reserve(_m.n_total_vertices() * 3);
    for(int vid = 0; vid < _m.n_total_vertices(); ++vid)
    {
        Mesh_connectivity::Vertex_iterator v = _m.vertex_at(vid);
        xyz.push_back(v.xyz().x());
        xyz.push_back(v.xyz().y());
        xyz.push_back(v.xyz().z());
    }

    _m.build_from_triangles(xyz, new_triangles);

    return true;
}


bool Loop_subdivider::subdivide_loop()
{
    using namespace mohecore;

    std::map<int, Mesh_connectivity::Vertex_iterator> edge_midpoints;

    // create new vertices at midpoints
    int n_old_vertices = _m.n_total_vertices();
    for (int heid = 0; heid < _m.n_total_half_edges(); ++heid)
    {
        auto he = _m.half_edge_at(heid);

        if (edge_midpoints.count(he.index()) > 0) continue;
        auto twin = he.twin();

        Eigen::Vector3d v0 = he.origin().xyz();
        Eigen::Vector3d v1 = he.dest().xyz();
        auto v2 = he.next().dest();
        auto v3 = twin.next().dest();
        
        Eigen::Vector3d new_pos =
            0.375 * (v0 + v1) + 0.125 * (v2.xyz() + v3.xyz());
            // Eigen::Vector3d new_pos =
            // 0.5 * (v0 + v1);

        // register new vertex
        auto new_vertex = _m.add_vertex();
        new_vertex.data().xyz = new_pos;

        // update half edge
        edge_midpoints[he.index()]   = new_vertex;
        edge_midpoints[twin.index()] = new_vertex;
    }

    // move old vertices
    std::vector<Eigen::Vector3d> updated_positions(n_old_vertices);
    for (int vid = 0; vid < n_old_vertices; ++vid)
    {
        auto v = _m.vertex_at(vid);
        if (!v.is_active()) continue;

        std::vector<Eigen::Vector3d> nbrs;
		auto ring = _m.vertex_ring_at(vid);
		do {
			auto neighbor = ring.half_edge().origin();
			if (neighbor.is_active() && neighbor.index() != vid) {
				nbrs.push_back(neighbor.xyz());
			}
		} while (ring.advance());

        int n = static_cast<int>(nbrs.size());
        // n = 6;
        double theta = 2.0 * M_PI / n;
        double w = (64.0 * n) / (40.0 - std::pow(3.0 + 2.0 * std::cos(2.0 * M_PI / n), 2.0)) - n;
        
        Eigen::Vector3d new_pos = Eigen::Vector3d::Zero();
        for (auto &p : nbrs) new_pos += p * 1 / (w + n);
        new_pos += v.xyz() * w / (w + n);

        updated_positions[vid] = new_pos;
    }

    // 
    std::vector<int> new_tris;
    for (int fid = 0; fid < _m.n_total_faces(); ++fid)
    {
        auto face = _m.face_at(fid);
        if (!face.is_active()) continue;

        auto he0 = face.half_edge();
        auto he1 = he0.next();
        auto he2 = he0.prev();

        int v0 = he0.origin().index();
        int v1 = he1.origin().index();
        int v2 = he2.origin().index();

        int m0 = edge_midpoints[he0.index()].index();
        int m1 = edge_midpoints[he1.index()].index();
        int m2 = edge_midpoints[he2.index()].index();

        new_tris.insert(new_tris.end(), { v0, m0, m2,
                                          v1, m1, m0,
                                          v2, m2, m1,
                                          m0, m1, m2 });
        face.deactivate();
    }

    // 4) reconstruct mesh
    std::vector<double> xyz;
    xyz.reserve(_m.n_total_vertices() * 3);

    for (int vid = 0; vid < _m.n_total_vertices(); ++vid)
    {
        auto v = _m.vertex_at(vid);
        if (!v.is_active()) continue;

        Eigen::Vector3d pos;
        if (vid < (int)updated_positions.size() && updated_positions[vid] != Eigen::Vector3d::Zero())
            pos = updated_positions[vid];
        else
            pos = v.xyz();

        xyz.push_back(pos.x());
        xyz.push_back(pos.y());
        xyz.push_back(pos.z());
    }

    _m.build_from_triangles(xyz, new_tris);

    return true;
}

} // end of mohecore
} // end of minimesh

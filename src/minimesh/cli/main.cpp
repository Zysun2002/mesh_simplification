//
// This is a bare executable that gets linked to the core library.
// You can use it to test your code, or start learning about the library.
//
// Here I have put an executable which reads a predefined mesh, flips a specific edge in that 
// mesh, and then writes it back to .vtk and .obj formats. Feel free to play with the example, change
// it, or move it to a completely different file.
//

#include <cstdio>
#include <iostream>

#define M_PI 3.14159265358979323846
#include <cmath>
#include <map>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;
fs::path out_path = "E:/Ziyu/workspace/course/geometryModeling/assignment_2/ZIyuSun_a2/mesh/output";
fs::path mesh_path = "E:/Ziyu/workspace/course/geometryModeling/assignment_2/ZIyuSun_a2/mesh";

#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/macros.hpp>

#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_simplification.hpp>

using namespace minimesh;

// ===================
// EXAMPLE UTILITY FUNCTIONS
// ===================

namespace 
{

//
// Create an example mesh file that we can read later.
//
void write_sphere_mesh()
{
    FILE *fl = fopen((out_path / "example_sphere.obj").string().c_str(), "w");
    if (!fl) return;

    // Vertices (icosahedron normalized)
    fprintf(fl, "v -0.525731 0.850651 0.000000\n");
    fprintf(fl, "v 0.525731 0.850651 0.000000\n");
    fprintf(fl, "v -0.525731 -0.850651 0.000000\n");
    fprintf(fl, "v 0.525731 -0.850651 0.000000\n");
    fprintf(fl, "v 0.000000 -0.525731 0.850651\n");
    fprintf(fl, "v 0.000000 0.525731 0.850651\n");
    fprintf(fl, "v 0.000000 -0.525731 -0.850651\n");
    fprintf(fl, "v 0.000000 0.525731 -0.850651\n");
    fprintf(fl, "v 0.850651 0.000000 -0.525731\n");
    fprintf(fl, "v 0.850651 0.000000 0.525731\n");
    fprintf(fl, "v -0.850651 0.000000 -0.525731\n");
    fprintf(fl, "v -0.850651 0.000000 0.525731\n\n");

    // Faces
    const int faces[20][3] = {
        {1,12,6}, {1,6,2}, {1,2,8}, {1,8,11}, {1,11,12},
        {2,6,10}, {6,12,5}, {12,11,3}, {11,8,7}, {8,2,9},
        {4,10,5}, {4,5,3}, {4,3,7}, {4,7,9}, {4,9,10},
        {5,10,6}, {3,5,12}, {7,3,11}, {9,7,8}, {10,9,2}
    };

    for (int i = 0; i < 20; ++i)
        fprintf(fl, "f %d %d %d\n", faces[i][0], faces[i][1], faces[i][2]);

    fclose(fl);
}

} // end of anonymus namespace

int main(int argc, char **argv)
{

  mohecore::Mesh_connectivity mesh;
  mohecore::Mesh_io io(mesh);
  mohecore::Simplifier simp(mesh);
  // mohecore::Loop_subdivider subdiv(mesh);


  printf("reading example mesh \n");
  io.read_obj_general((mesh_path / "cow_head.obj").string().c_str());


  int writing_index = 0;
  auto check_sanity_and_write_mesh = [&io, &mesh, &writing_index]()
  {
    force_assert( mesh.check_sanity_slowly() );

    printf("writing out_%d.obj \n", writing_index);

    io.write_obj( MINIMESH_STR((out_path / ("out_" + std::to_string(writing_index) + ".obj")).string()) );
    
    ++writing_index;
  };

// //   // Now check that the mesh is sane and write it  in both 
// //   // .vtk and .obj formats
//   // check_sanity_and_write_mesh();

// //   // Flip the edge between vertices 4 and 5 (the diagonal from lower right to upper left)
// //   // Note that the indices should become 0-index based rather than 1-index based.
// //   printf("subdivision \n");
//   // subdiv.subdivide_loop();
  simp.simplify_once();
  check_sanity_and_write_mesh();

  // simp.simplify();
  // check_sanity_and_write_mesh();

//   // subdiv.subdivide_loop();
//   // check_sanity_and_write_mesh();

//   // subdiv.subdivide_loop();
//   // check_sanity_and_write_mesh();


// //   modi.flip_edge( modi.get_halfedge_between_vertices(4-1, 5-1) );
// //   check_sanity_and_write_mesh();

// //   // Now flip back the edge again
// //   printf("flipping edge again ...\n");
// //   modi.flip_edge( modi.get_halfedge_between_vertices(2-1, 7-1) );
// //   check_sanity_and_write_mesh();

  return 0;
} // end of main()

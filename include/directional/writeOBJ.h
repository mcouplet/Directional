// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_WRITE_OBJ_H
#define DIRECTIONAL_WRITE_OBJ_H


#include <Eigen/Core>
#include <string>
#include <vector>
#include <directional/TriMesh.h>

namespace directional
{

/***Writing an OBJ file
Input:
fileName:     Name of file
mesh:         Surface mesh
TC:           Texture coordinates
FTC:          Per face (so |F|x3) indices into TC to associate texture to corners
mtlFileName:  name of MTL file that associates with the texture
textureName:  texture name inside the MTL file.
***/

bool writeOBJ(const std::string& fileName,
              const directional::TriMesh& mesh,
              const Eigen::MatrixXd& TC,
              const Eigen::MatrixXi& FTC,
              const std::string& mtlFileName = "",
              const std::string& textureName = "")
{
    std::ofstream out(fileName);
    if (!out.is_open())
        throw std::runtime_error("Failed to open file: " + fileName);
    
    // Write reference to MTL file if provided
    if (!mtlFileName.empty())
        out << "mtllib " << mtlFileName << "\n";
    
    // Use material if name is provided
    if (!textureName.empty())
        out << "usemtl " << textureName << "\n";
    
    // Write vertices
    for (int i = 0; i < mesh.V.rows(); ++i)
        out << "v " << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << "\n";
    
    // Write texture coordinates
    for (int i = 0; i < TC.rows(); ++i)
        out << "vt " << TC(i, 0) << " " << TC(i, 1) << "\n";
    
    // Write faces with texture indices
    for (int i = 0; i < mesh.F.rows(); ++i) {
        out << "f";
        for (int j = 0; j < 3; ++j) {
            int v_idx = mesh.F(i, j) + 1;    // OBJ uses 1-based indexing
            int vt_idx = FTC(i, j) + 1;
            out << " " << v_idx << "/" << vt_idx;
        }
        out << "\n";
    }
    
    out.close();
    return true;
}

bool inline writeOBJWithUV(
    const std::string& objFileName,
    const directional::TriMesh& mesh,
    const Eigen::MatrixXd& uv
) {

    int N = uv.cols() / 3;

    const Eigen::MatrixXd &V = mesh.V;
    const Eigen::MatrixXi &F = mesh.F;

    std::ofstream file(objFileName);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << objFileName << std::endl;
        return false;
    }

    // Sanity checks
    if (V.cols() != 3 || F.cols() != 3) {
        std::cerr << "Invalid V or F dimensions." << std::endl;
        return false;
    }
    if (uv.rows() != F.rows()) {
        std::cerr << "uv must have one row per face." << std::endl;
        return false;
    }
    if (uv.cols() < 6) {
        std::cerr << "uv must have at least 6 columns (3 corners × 2)." << std::endl;
        return false;
    }

    // Write vertices
    for (int i = 0; i < V.rows(); ++i) {
        file << "v "
             << V(i,0) << " "
             << V(i,1) << " "
             << V(i,2) << "\n";
    }

    // Write texture coordinates (one per face corner)
    // vt indices will be: 1, 2, 3, ..., 3*F.rows()
    for (int f = 0; f < F.rows(); ++f) {
        for (int c = 0; c < 3; ++c) {
            double u = uv(f, N*c + 0);
            double v = uv(f, N*c + 1);
            file << "vt " << u << " " << v << "\n";
        }
    }

    // Write faces
    for (int f = 0; f < F.rows(); ++f) {
        file << "f ";
        for (int c = 0; c < 3; ++c) {
            int v_idx  = F(f, c) + 1;          // OBJ is 1-indexed
            int vt_idx = 3*f + c + 1;          // per-corner vt
            file << v_idx << "/" << vt_idx;
            if (c < 2) file << " ";
        }
        file << "\n";
    }

    // Write feature edges as line primitives
    if (mesh.isFeatureEdge.size() == mesh.EV.rows()) {
        for (int e = 0; e < mesh.EV.rows(); ++e) {
            if (mesh.isFeatureEdge[e]) {
                int v1 = mesh.EV(e, 0) + 1;  // OBJ is 1-indexed
                int v2 = mesh.EV(e, 1) + 1;
                file << "l " << v1 << " " << v2 << "\n";
            }
        }
    }

    file.close();
    return true;
}


// writes a polygonal mesh as an ascii OBJ file
// Inputs:
//   str  path to .obj file
//   V    eigen double matrix  #V by 3 - vertex coordinates
//   D    eigen int vector     #F by 1 - face degrees
//   F    eigen int matrix     #F by max(D) - vertex indices in face
inline bool polygonal_write_OBJ(const std::string& str,
                                const Eigen::MatrixXd& V,
                                const Eigen::VectorXi& D,
                                const Eigen::MatrixXi& F)
{
  using namespace std;
  using namespace Eigen;

  ofstream file(str);
  if (!file.is_open())
    return false;

  file.precision(9);
  file << std::fixed;

  // Write vertices
  for (int i = 0; i < V.rows(); ++i)
  {
    file << "v "
         << V(i, 0) << " "
         << V(i, 1) << " "
         << V(i, 2) << "\n";
  }

  // Write faces (OBJ is 1-based indexing)
  for (int i = 0; i < F.rows(); ++i)
  {
    file << "f";
    for (int j = 0; j < D(i); ++j)
    {
      file << " " << (F(i, j) + 1);
    }
    file << "\n";
  }

  file.close();
  return true;
}

}

#endif

// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_READ_OBJ_H
#define DIRECTIONAL_READ_OBJ_H
#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <directional/TriMesh.h>


namespace directional
{

//Reading an OBJ file into a surface TriMesh class.
bool inline readOBJ(const std::string objFileName,
                    directional::TriMesh& mesh) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::ifstream file(objFileName);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << objFileName << std::endl;
        return false;
    }
    
    std::string line;
    std::vector<Eigen::RowVector3d> vertexList;
    std::vector<Eigen::RowVector3i> faceList;
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string type;
        ss >> type;
        
        if (type == "v") {
            Eigen::RowVector3d newVertex;
            ss >> newVertex(0) >> newVertex(1) >> newVertex(2);
            vertexList.push_back(newVertex);
        } else if (type == "f") {
            Eigen::RowVector3i newFace;
            std::string v_str;
            
            for (int i = 0; i < 3; ++i) {
                ss >> v_str;
                std::stringstream vs(v_str);
                std::string idx_str;
                std::getline(vs, idx_str, '/'); // take only the first part
                newFace(i) = std::stoi(idx_str); // OBJ is 1-indexed
            }
            
            faceList.push_back(newFace);
        }
    }
    
    file.close();
    
    V.resize(vertexList.size(), 3);
    for (int i = 0; i < vertexList.size(); ++i)
        V.row(i) = vertexList[i];
    
    F.resize(faceList.size(), 3);
    for (int i = 0; i < faceList.size(); ++i)
        F.row(i) = faceList[i];
    
    // Convert to 0-indexing
    int minIndex = F.minCoeff();
    F.array() -= minIndex;
    
    mesh.set_mesh(V, F);
    return true;
}

bool inline readOBJWithUV(
    const std::string& objFileName,
    directional::TriMesh& mesh,
    Eigen::MatrixXd& uv
) {
    std::ifstream file(objFileName);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << objFileName << std::endl;
        return false;
    }

    std::string line;
    std::vector<Eigen::RowVector3d> vertexList;
    std::vector<Eigen::RowVector2d> uvList;
    std::vector<Eigen::RowVector3i> faceList;
    std::vector<Eigen::RowVector3i> faceUVList;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string type;
        ss >> type;

        if (type == "v") {
            Eigen::RowVector3d v;
            ss >> v(0) >> v(1) >> v(2);
            vertexList.push_back(v);
        }
        else if (type == "vt") {
            Eigen::RowVector2d t;
            ss >> t(0) >> t(1);
            uvList.push_back(t);
        }
        else if (type == "f") {
            Eigen::RowVector3i f, ft;
            std::string token;

            for (int i = 0; i < 3; ++i) {
                ss >> token;

                std::stringstream ts(token);
                std::string v_str, vt_str;

                std::getline(ts, v_str, '/');
                std::getline(ts, vt_str, '/'); // may be empty if no vt

                f(i) = std::stoi(v_str);
                ft(i) = vt_str.empty() ? -1 : std::stoi(vt_str);
            }

            faceList.push_back(f);
            faceUVList.push_back(ft);
        }
    }

    file.close();

    // Build V
    Eigen::MatrixXd V;
    V.resize(vertexList.size(), 3);
    for (int i = 0; i < (int)vertexList.size(); ++i)
        V.row(i) = vertexList[i];

    // Build F
    Eigen::MatrixXi F;
    F.resize(faceList.size(), 3);
    for (int i = 0; i < (int)faceList.size(); ++i)
        F.row(i) = faceList[i];

    // Convert to 0-indexing
    int minV = F.minCoeff();
    F.array() -= minV;

    // Build uv as F x 6 (per-face, per-corner)
    uv.resize(faceUVList.size(), 6);
    uv.setZero();

    for (int f = 0; f < (int)faceUVList.size(); ++f) {
        for (int c = 0; c < 3; ++c) {
            int tid = faceUVList[f](c);
            if (tid >= 0) {
                const auto& t = uvList[tid - 1]; // OBJ is 1-indexed
                uv(f, 2*c + 0) = t(0);
                uv(f, 2*c + 1) = t(1);
            }
        }
    }

    mesh.set_mesh(V, F);
    return true;
}
}

#endif

// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can

// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CURL_MATCHING_H
#define DIRECTIONAL_CURL_MATCHING_H

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/effort_to_indices.h>

namespace directional
{
// Takes a field in raw form and computes both the curl-matching effort and the consequent curl matching on every tangent-space adjacency.
// Important: if the Raw field in not CCW ordered, the result is meaningless.
// Note: curl is only (future work...) defined for face-based fields.
// Input:
//  rawField:   RAW_FIELD type field
// Output:
//  curlNorm:   L2-norm of the curl vector, per edge
//  rawField:   With input field matching, effort, and singularities computed
inline void curl_matching(directional::CartesianField& rawField,
                          Eigen::VectorXd& curlNorm)
{
    
    typedef std::complex<double> Complex;
    using namespace Eigen;
    using namespace std;
    
    //this only works on face-based fields for now
    assert(rawField.tb->discTangType()==discTangTypeEnum::FACE_SPACES && "This function only supports face-based fields for now.");
    PCFaceTangentBundle* ftb = (PCFaceTangentBundle*)(rawField.tb);
    rawField.matching.conservativeResize(ftb->mesh->EF.rows());
    rawField.matching.setConstant(-1);
    curlNorm = Eigen::VectorXd::Zero(ftb->mesh->EF.rows());
    
    // Compute edge vectors, i.e., vector = vertex2 - vertex1, normalized
    MatrixXd edgeVectors(ftb->mesh->EF.rows(), 3);
    for (int i = 0; i < ftb->mesh->EF.rows(); i++) {
        if (ftb->mesh->EF(i, 0) == -1 || ftb->mesh->EF(i, 1) == -1)
            continue;
        edgeVectors.row(i) = (ftb->mesh->V.row(ftb->mesh->EV(i, 1)) - ftb->mesh->V.row(ftb->mesh->EV(i, 0))).normalized();
    }
    
    //effort = VectorXd::Zero(rawField.mesh->EF.rows());
    for (int i = 0; i < ftb->mesh->EF.rows(); i++) { // for every edge
        if (ftb->mesh->EF(i, 0) == -1 || ftb->mesh->EF(i, 1) == -1)
            continue;
        //computing free coefficient effort (a.k.a. [Diamanti et al. 2014])
        //Complex freeCoeffEffort(1.0, 0.0);
        int indexMinFromZero=0;
        //finding where the 0 vector in EF(i,0) goes to with smallest rotation angle in EF(i,1), computing the effort, and then adjusting the matching to have principal effort.
        // Matteo's note: This is curl-based matching! It's not the same as smoothness-based matching.
        // So the topology might be different from the singularity structure which uses smoothness.
        double minCurl = 32767000.0;
        for (int j = 0; j < rawField.N; j++) { // j is the candidate matching
            double currCurl = 0;
            for (int k=0;k<rawField.N;k++){
                RowVector3d vecDiff =rawField.extField.block(ftb->mesh->EF(i, 1), 3 * ((j+k)%rawField.N), 1, 3)-rawField.extField.block(ftb->mesh->EF(i, 0), 3*k, 1, 3); // subtract extrinsic vectors
                // Shouldn't the triangles be unfolded here??
                currCurl +=pow(edgeVectors.row(i).dot(vecDiff),2.0);
            }
            
            if (currCurl < minCurl){
                indexMinFromZero=j;
                minCurl=currCurl;
            }
        }
        
        // The matching is the one that minimizes the total squared curl of matched vector fields
        rawField.matching(i) =indexMinFromZero;
        curlNorm(i)= sqrt(minCurl);
        
        //computing the full effort for 0->indexMinFromZero, and readjusting the matching to fit principal effort
        // where is this readjusting done??
        // Matteo's note: effort and freeCoeff do NOT depend on the matching. Need to understand that better
        Complex freeCoeff(1,0);
        rawField.effort.resize(rawField.matching.size());
        for (int j = 0; j < rawField.N; j++) { // for every frame vector j
            // Two neighboring faces: f and g
            //RowVector3d vecjf = rawField.extField.block(rawField.adjSpaces(i, 0), 3*j, 1, 3);
            RowVector2d vecjf = rawField.intField.block(ftb->adjSpaces(i, 0), 2 * j, 1, 2);
            Complex vecjfc = Complex(vecjf(0),vecjf(1));
            RowVector2d vecjg = rawField.intField.block(ftb->adjSpaces(i, 1), 2 * ((rawField.matching(i)+j+rawField.N)%rawField.N), 1, 2);
            Complex vecjgc = Complex(vecjg(0),vecjg(1));
            Complex transvecjfc = vecjfc*ftb->connection(i); // vecjfc transported to face g
            freeCoeff *= (vecjgc / transvecjfc); // the effort. Note that arg(u*v) = arg(u) + arg(v) (mod 2π)

        }
        
        rawField.effort(i) = arg(freeCoeff);

        // Matteo's version: pick the matching that minimizes the sum of angles.
        // This will overwrite rawField's matching
        // I think it's a way to do smoothness-based matching
        double minSumArgs = DBL_MAX, bestMatching = -1;
        for (int m = 0; m < rawField.N; m++) {
            Complex freeCoeff(1,0);
            double sumArgs = 0;
            for (int j = 0; j < rawField.N; j++) { // for every frame vector j
                // Two neighboring faces: f and g
                RowVector2d vecjf = rawField.intField.block(ftb->adjSpaces(i, 0), 2 * j, 1, 2);
                Complex vecjfc = Complex(vecjf(0),vecjf(1));
                RowVector2d vecjg = rawField.intField.block(ftb->adjSpaces(i, 1), 2 * ((m+j+rawField.N)%rawField.N), 1, 2);
                Complex vecjgc = Complex(vecjg(0),vecjg(1));
                Complex transvecjfc = vecjfc*ftb->connection(i); // vecjfc transported to face g
                freeCoeff *= (vecjgc / transvecjfc); // the effort
                sumArgs += abs(arg(vecjgc / transvecjfc)); // this can be seen as a smoothness measure
            }

            if (sumArgs < minSumArgs) {
                bestMatching = m;
                minSumArgs = sumArgs;
                rawField.matching(i) = bestMatching;
            }
        }
    }

    Eigen::VectorXd indices = rawField.tb->cycles * rawField.matching.cast<double>();

    //Getting final singularities and their indices
    effort_to_indices(rawField);
    
}
}


#endif



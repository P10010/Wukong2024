#include "../include/PBD.h"
#include <cmath>
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/triangle_triangle_adjacency.h>

#include "../include/Util.h"
void PBD::initializeFromFile(const std::string& filename)
{
    scene = filename.data();
    MatrixXT V;
    MatrixXi F;
    igl::readOBJ(scene, V, F);

    // TODO: necessary?
    // these functions copy the V, F, that is, what was read on the obj, to the
    // internal representation atRest and faces they were copied from
    // DiscreteShell and also "fatten" it, but that is probably unnecessary for
    // us. you could set the second parameter to 1 rather than 3 if you didn't
    // want to fatten (but at that point it might be worth it creating another
    // function in Util.h) you would probably need to change app.initializeScene
    // too maybe fatten is actually necessary to give some "thickness" to the
    // surface mesh? idk. Changing it breaks something else, I'll leave it as is
    // for now
    // iglMatrixFatten<T, 3>(V, atRest);
    // iglMatrixFatten<int, 3>(F, faces);

    faces = F;
    atRest = V;
    currentV = atRest;
    int nRows = currentV.rows();
    p.resize(nRows, 3);
    w.resize(nRows);
    // TODO: put 1/mass, at the moment I set the mass to be 1 for every vertex
    // for every vertex i, it should be 1/3 of the mass of each adjacent
    // triangle
    for (int i = 0; i < nRows; i++)
    {
        w(i) = 1;
    }

    igl::edges(F, edges);
    edge_lengths.resize(edges.rows());
    for (int i = 0; i < edges.rows(); i++)
    {
        int u = edges(i, 0);
        int v = edges(i, 1);
        edge_lengths(i) = (atRest.row(u) - atRest.row(v)).norm();
    }

    MatrixXT normals(faces.rows(), 3);
    igl::per_face_normals(atRest, faces, normals);
    normals.rowwise().normalize();

    igl::triangle_triangle_adjacency(F, TT, TTi);

    dihedral_angles.resize(faces.rows(), 3);
    for (int i = 0; i < faces.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (TT(i, j) < 0)
                continue;
            TV n1, n2;
            n1 = normals.row(i);
            n2 = normals.row(TT(i, j));
            dihedral_angles(i, j) = std::acos(n1.dot(n2));
        }
    }

    v = MatrixXT::Zero(nRows, 3);

    // pinned diagonal
    // posConstraintsIdxs.resize(2);
    // posConstraintsIdxs << 9, 90;
    // posConstraintsV.resize(2, 3);
    // posConstraintsV << -0.3, 0.0, 0.0, 0.0, 0.3, 0.0;

    // pinned horizontal
    posConstraintsIdxs.resize(2);
    posConstraintsIdxs << 0, 90;
    posConstraintsV.resize(2, 3);
    posConstraintsV << 0.0, 0.0, 0.0, 0.0, 0.3, 0.0;
}

void PBD::stretchingConstraints(int solver_it)
{
    for (int i = 0; i < edges.rows(); i++)
    {
        int a, b;
        a = edges(i, 0);
        b = edges(i, 1);
        T l0 = edge_lengths(i);

        TV p1, p2;
        p1 = p.row(a);
        p2 = p.row(b);

        TV dp1 = -w(a) / (w(a) + w(b)) * ((p1 - p2).norm() - l0) * (p1 - p2).normalized();
        TV dp2 = w(b) / (w(a) + w(b)) * ((p1 - p2).norm() - l0) * (p1 - p2).normalized();

        T alpha = 1. - std::pow(1. - k_stretch, 1. / solver_it);
        p.row(a) += dp1 * alpha;
        p.row(b) += dp2 * alpha;
    }
}

void PBD::bendingConstraints(int solver_it)
{
    for (int f = 0; f < faces.rows(); f++)
    {
        for (int e = 0; e < 3; e++)
        {
            int fp = TT(f, e); // adjacent face to fp opposite to edge e
            // skip adjacent face if we have already processed it (or if -1)
            if (fp < f) continue;

            int p1i, p2i, p3i, p4i;
            p1i = faces(f, (e + 2) % 3);
            p2i = faces(f, (e + 0) % 3);
            p3i = faces(f, (e + 1) % 3);
            for (int j = 0; j < 3; j++)
            {
                if (faces(fp, j) != p1i && faces(fp, j) != p2i && faces(fp, j) != p3i)
                {
                    p4i = faces(fp, j);
                    break;
                }
            }

            TV p1, p2, p3, p4;
            p1 = p.row(p1i);
            p2 = p.row(p2i);
            p3 = p.row(p3i);
            p4 = p.row(p4i);
            // substract p1 to make calculations easier
            p2 -= p1;
            p3 -= p1;
            p4 -= p1;

            TV n1, n2;
            n1 = p2.cross(p3).normalized();
            n2 = p2.cross(p4).normalized();

            T d = std::min(1.0, n1.dot(n2));

            TV q1, q2, q3, q4;
            q3 = (p2.cross(n2) + n1.cross(p2) * d) / p2.cross(p3).norm();
            q4 = (p2.cross(n1) + n2.cross(p2) * d) / p2.cross(p4).norm();
            q2 = -(p3.cross(n2) + n1.cross(p3) * d) / p2.cross(p3).norm() -
                 (p4.cross(n1) + n2.cross(p4) * d) / p2.cross(p4).norm();
            q1 = -q2 - q3 - q4;

            T sum_qnorm2 = q1.squaredNorm() + q2.squaredNorm() +
                           q3.squaredNorm() + q4.squaredNorm();
            if (sum_qnorm2 == 0) continue;
            T sum_w = w(p1i) + w(p2i) + w(p3i) + w(p4i);
            T sqrt_1md2 = std::sqrt(1.0 - d * d);
            T acosd = std::acos(d);
            T phi0 = dihedral_angles(f, e);

            T alpha = 1. - std::pow(1. - k_bend, 1. / solver_it);

            TV dp1, dp2, dp3, dp4;
            dp1 = -4 * w(p1i) / sum_w * sqrt_1md2 * (acosd - phi0) / sum_qnorm2 * q1;
            dp2 = -4 * w(p2i) / sum_w * sqrt_1md2 * (acosd - phi0) / sum_qnorm2 * q2;
            dp3 = -4 * w(p3i) / sum_w * sqrt_1md2 * (acosd - phi0) / sum_qnorm2 * q3;
            dp4 = -4 * w(p4i) / sum_w * sqrt_1md2 * (acosd - phi0) / sum_qnorm2 * q4;

            p.row(p1i) += dp1 * alpha;
            p.row(p2i) += dp2 * alpha;
            p.row(p3i) += dp3 * alpha;
            p.row(p4i) += dp4 * alpha;
        }
    }
}

void PBD::positionConstraints()
{
    for(int i = 0; i < posConstraintsIdxs.rows(); ++i)
    {
        int idx = posConstraintsIdxs[i];
        p.row(idx) = posConstraintsV.row(i);
    }
    for (int i = 0; i < p.rows(); ++i)
    {
        // floor at z=-1/4
        p(i, 2) = std::max(p(i, 2), -0.25);
    }
}

void PBD::projectConstraints(int solver_it)
{
    if (stretchingConstraintsActivated)    
        PBD::stretchingConstraints(solver_it);

    if (bendingConstraintsActivated)
        PBD::bendingConstraints(solver_it);

    // TODO: collision constraints

    if (positionConstraintsActivated)
        PBD::positionConstraints();
}

bool PBD::advanceOneStep(int step)
{

    // Pre-solve: apply external forces
    int nRows = currentV.rows();
    for (int i = 0; i < nRows; i++)
    {
        v.row(i) += dt * g;
    }

    // TODO: damp velocities? potentially not necessary?

    for (int i = 0; i < nRows; i++)
    {
        // Verlet integration, as in paper
        p.row(i) = currentV.row(i) + dt * v.row(i);
    }

    // generate collision constraints

    // Solve constraints
    for (int it = 0; it < numIterations; it++)
    {
        projectConstraints(it);
    }

    // Post-solve: update velocities
    for (int i = 0; i < nRows; i++)
    {
        v.row(i) << (p.row(i) - currentV.row(i)) * 1. / dt;
        currentV.row(i) = p.row(i);
    }

    // TODO: velocity update


    if (step == 1000)
        return true;
    return false;
}

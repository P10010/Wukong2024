#include "../include/PBD.h"
#include <cmath>
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/triangle_triangle_adjacency.h>
#include <random>

void PBD::initializeFromFile(const std::string& filename)
{
    scene = filename.data();
    MatrixXT V;
    MatrixXi F;
    igl::readOBJ(scene, V, F);

    // VectorXT offset = VectorXT::Zero(V.cols());
    // offset[0] = 0.0; 
    // offset[1] = 3.13 - 3.22;
    // offset[2] = 4.56 - 4.5;
    //
    // std::cout << "pre  V.row(59) = " << V.row(59)<< std::endl;
    // for (int i = 0; i < V.rows(); ++i) {
    //     V(i,0) += 10;
    //     V(i,1) += 10;
    //     V(i,2) += 10;
    // }
    // std::cout << "post V(59,0) = " << V.row(59)<< std::endl;

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

    // normalize the mesh to fit in a unit cube
    // for (int i = 0; i < V.rows(); i++)
    // {
    //     V.row(i) -= V.colwise().minCoeff();
    //     V.row(i) /= V.colwise().maxCoeff().maxCoeff();
    // }

    faces = F;
    atRest = V;
    currentV = atRest;
    int nRows = currentV.rows();
    p.resize(nRows, 3);
    m.resize(nRows);
    w.resize(nRows);

    // Calculate mass for each vertex based on adjacent triangles
    for (int i = 0; i < nRows; i++)
    {
        m(i) = 0;
    }
    for (int i = 0; i < faces.rows(); i++)
    {
        int v1 = faces(i, 0);
        int v2 = faces(i, 1);
        int v3 = faces(i, 2);
        TV a = atRest.row(v2) - atRest.row(v1);
        TV b = atRest.row(v3) - atRest.row(v1);
        T area = a.cross(b).norm() * 0.5;
        m(v1) += area / 3.0;
        m(v2) += area / 3.0;
        m(v3) += area / 3.0;
    }

    for (int i = 0; i < nRows; i++)
    {
        w(i) = 1. / m(i);
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

    // pinned to mast (sail)
    posConstraintsIdxs = -1 * VectorXi::Ones(6);
    posConstraintsV.resize(6, 3);
    posConstraintsIdxs[0] = 1;
    //
    posConstraintsIdxs[1] = 50;
    posConstraintsIdxs[2] = 57;
    //
    posConstraintsIdxs[3] = 3;
    posConstraintsIdxs[4] = 0;
    posConstraintsIdxs[5] = 2;
    for (int i = 0; i < posConstraintsIdxs.rows(); ++i)
    {
        int idx = posConstraintsIdxs[i];
        posConstraintsV.row(i) = atRest.row(idx);
        w(idx) = 0;
    }

    lambdas = VectorXT::Zero(constraint_idx.size());

    dist10 = std::uniform_int_distribution<int>(0, 10);
    distV = std::uniform_int_distribution<int>(0, V.rows() - 1);
    distF = std::uniform_int_distribution<int>(0, F.rows() - 1);
}

void PBD::stretchingConstraintsXPBD()
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

        T w_sum = w(a) + w(b);

        T alpha = alpha_stretch / (dt * dt);
        T lambda = lambdas[constraint_idx["stretch"]];

        if (w_sum == 0) continue;
        // XXX:
        // T dlambda = -((p1 - p2).norm() - l0 + alpha * lambda) / (w_sum + alpha);
        T dlambda = -((p1 - p2).norm() - l0) / (w_sum + alpha);

        // if (dlambda > 1.0) {
        //     std::cout << "===================" << std::endl;
        //     std::cout << w_sum << std::endl;
        //     std::cout << alpha << std::endl;
        //     std::cout << lambda << std::endl;
        // }

        TV dp1 = w(a) * dlambda * (p1 - p2).normalized();
        TV dp2 = w(b) * dlambda * (p2 - p1).normalized();

        // if (dp1.norm() > 1.0 or dp2.norm() > 1.0) {
        //     std::cout << "--------------------" << std::endl;
        //     std::cout << "a: " << a << std::endl;
        //     std::cout << "|dp1(a)|: " << dp1.norm() << std::endl;
        //     std::cout << "b: " << b << std::endl;
        //     std::cout << "|dp2(b)|: " << dp2.norm() << std::endl;
        //     std::cout << "w_sum: " << w_sum << std::endl;
        //     std::cout << "alpha: " << alpha << std::endl;
        //     std::cout << "lambda: " << lambda << std::endl;
        // }

        p.row(a) += dp1;
        p.row(b) += dp2;

        lambdas[constraint_idx["stretch"]] += dlambda;
    }
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

        TV dp1 = -w(a) / (w(a) + w(b)) * ((p1 - p2).norm() - l0) *
                 (p1 - p2).normalized();
        TV dp2 = w(b) / (w(a) + w(b)) * ((p1 - p2).norm() - l0) *
                 (p1 - p2).normalized();

        T alpha = 1. - std::pow(1. - k_stretch, 1. / solver_it);
        p.row(a) += dp1 * alpha;
        p.row(b) += dp2 * alpha;
    }
}

void PBD::bendingConstraintsXPBD()
{
    for (int f = 0; f < faces.rows(); f++)
    {
        for (int e = 0; e < 3; e++)
        {
            int fp = TT(f, e); // adjacent face to fp opposite to edge e
            // skip adjacent face if we have already processed it (or if -1)
            if (fp < f)
                continue;

            int p1i, p2i, p3i, p4i;
            p1i = faces(f, (e + 2) % 3);
            p2i = faces(f, (e + 0) % 3);
            p3i = faces(f, (e + 1) % 3);
            for (int j = 0; j < 3; j++)
            {
                if (faces(fp, j) != p1i && faces(fp, j) != p2i &&
                    faces(fp, j) != p3i)
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

            T sum_wqnorm2 = w(p1i) * q1.squaredNorm() + w(p2i) * q2.squaredNorm() +
                           w(p3i) * q3.squaredNorm() + w(p4i) * q4.squaredNorm();
            T acosd = std::acos(d);
            T phi0 = dihedral_angles(f, e);

            T alpha = alpha_bend / (dt * dt);
            T lambda = lambdas[constraint_idx["bend"]];

            if (sum_wqnorm2 == 0) continue;
            // XXX:
            // T dlambda = -(acosd - phi0 + alpha * lambda) / (sum_wqnorm2 + alpha);
            T dlambda = -(acosd - phi0) / (sum_wqnorm2 + alpha);

            T sqrt_1md2 = std::sqrt(1.0 - d * d);

            TV dp1, dp2, dp3, dp4;
            dp1 = w(p1i) * dlambda * sqrt_1md2 * q1;
            dp2 = w(p2i) * dlambda * sqrt_1md2 * q2;
            dp3 = w(p3i) * dlambda * sqrt_1md2 * q3;
            dp4 = w(p4i) * dlambda * sqrt_1md2 * q4;

            p.row(p1i) += dp1;
            p.row(p2i) += dp2;
            p.row(p3i) += dp3;
            p.row(p4i) += dp4;

            lambdas[constraint_idx["bend"]] += dlambda;
        }
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
            if (fp < f)
                continue;

            int p1i, p2i, p3i, p4i;
            p1i = faces(f, (e + 2) % 3);
            p2i = faces(f, (e + 0) % 3);
            p3i = faces(f, (e + 1) % 3);
            for (int j = 0; j < 3; j++)
            {
                if (faces(fp, j) != p1i && faces(fp, j) != p2i &&
                    faces(fp, j) != p3i)
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
            if (sum_qnorm2 == 0)
                continue;
            T sum_w = w(p1i) + w(p2i) + w(p3i) + w(p4i);
            T sqrt_1md2 = std::sqrt(1.0 - d * d);
            T acosd = std::acos(d);
            T phi0 = dihedral_angles(f, e);

            T alpha = 1. - std::pow(1. - k_bend, 1. / solver_it);

            TV dp1, dp2, dp3, dp4;
            dp1 = -4 * w(p1i) / sum_w * sqrt_1md2 * (acosd - phi0) /
                  sum_qnorm2 * q1;
            dp2 = -4 * w(p2i) / sum_w * sqrt_1md2 * (acosd - phi0) /
                  sum_qnorm2 * q2;
            dp3 = -4 * w(p3i) / sum_w * sqrt_1md2 * (acosd - phi0) /
                  sum_qnorm2 * q3;
            dp4 = -4 * w(p4i) / sum_w * sqrt_1md2 * (acosd - phi0) /
                  sum_qnorm2 * q4;

            p.row(p1i) += dp1 * alpha;
            p.row(p2i) += dp2 * alpha;
            p.row(p3i) += dp3 * alpha;
            p.row(p4i) += dp4 * alpha;
        }
    }
}

// check whether point q intersects the triangle with vertices p1, p2, p3
bool PBD::pointIntersectsTriangle(const PBD::TV& q, const PBD::TV& p1,
                                  const PBD::TV& p2, const PBD::TV& p3) const
{
    // Compute triangle normal
    TV cr = (p2 - p1).cross(p3 - p1);
    TV n = cr.normalized();
    T dist = (q - p1).dot(n);

    int fromBelow;
    // check if vertex is from below with respect to the triangle normal (it
    // should be correct)
    if (dist < 0)
    {
        fromBelow = -1;
    }
    else
        fromBelow = 1;

    // note: if dist is small the point is close to the plane containing the
    // triangle, but I need to do some extra checks to verify that is actually
    // close to the triangle
    if (fromBelow * dist > h)
        return false;

    // calculate barycentric coordinates of projection into the plane (from
    // Robust Treatment of Collisions, Contact and Friction for Cloth Animation)
    TV x13 = p1 - p3, x23 = p2 - p3, x43 = q - p3;
    TM2 A;
    A << x13.dot(x13), x13.dot(x23), x13.dot(x23), x23.dot(x23);
    TV2 b = {x13.dot(x43), x23.dot(x43)};
    TV2 w12 = A.colPivHouseholderQr().solve(b);
    T w3 = 1 - w12[0] - w12[1];
    // δ is h divided by a characteristic length of the triangle, i.e. squared
    // root of the area
    T area = 0.5 * cr.norm();
    if (area == 0)
        return false;
    T delta = h / (sqrt(abs(area)));
    return (w12.x() >= -delta && w12.x() <= 1 + delta && w12.y() >= -delta &&
            w12.y() <= 1 + delta && w3 >= -delta &&
            w3 <= 1 + delta); // collision detected
}

int PBD::hash(int i, int j, int k, int n)
{
    return (((long)i * prime1) ^ ((long)j * prime2) ^ ((long)k * prime3)) % n;
}
int PBD::hash(PBD::TV p, T l, int n, TV& minCoord)
{
    int i = std::floor((p.x() - minCoord.x()) / l),
        j = std::floor((p.y() - minCoord.y()) / l),
        k = std::floor((p.z() - minCoord.z()) / l);
    return hash(i, j, k, n);
}

void PBD::hashVertices(std::vector<std::vector<int>>& hashTable, T boxSize,
                       TV& minCoord)
{
    int n = hashTable.size();
    for (int i = 0; i < p.rows(); ++i)
    {
        TV q = p.row(i);
        int hashVal = hash(q, boxSize, n, minCoord);
        hashTable[hashVal].push_back(i);
        //    if(hashTable[hashVal].size()>10)
        //      std::cout<<hashTable[hashVal].size()<<"\n";
    }
}

void PBD::spatialHashing()
{
    TV minCoord = p.colwise().minCoeff();
    TV maxCoord = p.colwise().maxCoeff();

    // T coeff = 1; //(pow(p.rows(),T(0.33)));
    T boxSize = ((maxCoord - minCoord).maxCoeff()) / 100;
    // T boxSize = 2 * h;
    int n = 4 * p.rows();
    std::vector<std::vector<int>> hashTable(n);
    hashVertices(hashTable, boxSize, minCoord);

    for (int f = 0; f < faces.rows(); f++)
    {
        int t1 = faces(f, 0), t2 = faces(f, 1), t3 = faces(f, 2);
        TM triangleVertices;
        TV p1 = p.row(t1), p2 = p.row(t2), p3 = p.row(t3);
        triangleVertices << p1, p2, p3;
        std::array<int, 3> minBox, maxBox;
        // std::cout << "p1\n" << p1 << std::endl;
        // std::cout << "triangle vertices\n"
        //           << triangleVertices << std::endl;
        for (int i = 0; i < 3; i++)
        {
            TV tmp = triangleVertices.row(i);
            minBox[i] = (int)std::floor(
                (triangleVertices.row(i).minCoeff() - minCoord(i)) / boxSize);
            // std::max(0,(int)std::floor((triangleVertices.row(i).minCoeff()-minCoord(i))/boxSize)-1);
            maxBox[i] = (int)std::floor(
                (triangleVertices.row(i).maxCoeff() - minCoord(i)) / boxSize);
        }
        // std::cout << "minBox\n" << minBox[0] << " " << minBox[1] << " "
        //           << minBox[2] << std::endl;
        // std::cout << "maxBox\n" << maxBox[0] << " " << maxBox[1] << " "
        //             << maxBox[2] << std::endl;

        for (int i = minBox[0]; i <= maxBox[0]; i++)
            for (int j = minBox[1]; j <= maxBox[1]; j++)
                for (int k = minBox[2]; k <= maxBox[2]; k++)
                {
                    int hashVal = hash(i, j, k, n);
                    for (int qIdx : hashTable[hashVal])
                        if (qIdx != t1 && qIdx != t2 && qIdx != t3 &&
                            pointIntersectsTriangle(p.row(qIdx), p1, p2, p3))
                        {
                            CollisionConstraint constraint;
                            constraint.qIdx = qIdx;
                            constraint.f = faces.row(f);
                            collisionConstraintsList.push_back(constraint);
                        }
                }
    }
}
void PBD::generateCollisionConstraints()
{
    // Clear previous collision constraints
    collisionConstraintsList.clear();

    if (useSpatialHashing)
    {
        spatialHashing();
    }
    else
    {
        // Iterate through all vertices
        for (int i = 0; i < p.rows(); ++i)
        {
            // Iterate through all triangles
            for (int f = 0; f < faces.rows(); f++)
            {
                int t1 = faces(f, 0), t2 = faces(f, 1), t3 = faces(f, 2);
                if (t1 == i || t2 == i || t3 == i)
                    continue;
                TV q = p.row(i);
                TV p1 = p.row(t1), p2 = p.row(t2), p3 = p.row(t3);

                if (pointIntersectsTriangle(q, p1, p2, p3))
                {
                    CollisionConstraint constraint;
                    constraint.qIdx = i;
                    constraint.f = faces.row(f);
                    collisionConstraintsList.push_back(constraint);
                }
            }
        }
    }
}

void PBD::collisionConstraints()
{
    for (auto& constraint : collisionConstraintsList)
    {
        int t1 = constraint.f(0), t2 = constraint.f(1), t3 = constraint.f(2);

        TV q = p.row(constraint.qIdx);
        TV p1 = p.row(t1), p2 = p.row(t2), p3 = p.row(t3);

        // Project the vertex and triangle vertices to satisfy the constraint
        T w_q = w(constraint.qIdx);
        T w_p1 = w(t1), w_p2 = w(t2), w_p3 = w(t3);

        T sum_w = w_q + w_p1 + w_p2 + w_p3;
        if (sum_w == 0)
            continue;

        // Compute triangle normal
        TV cr = (p2 - p1).cross(p3 - p1);
        TV n = cr.normalized();
        T dist = (q - p1).dot(n);

        int fromBelow;
        // check if vertex is from below with respect to the triangle normal (it
        // should be correct)
        if (dist < 0)
            fromBelow = -1;
        else
            fromBelow = 1;

        if (fromBelow * dist < h)
        {
            // I manually calculated the derivative, I know this code is pretty
            // ugly, and probably it could be written in a more compact form.
            // Hopefully the calculations are correct
            TV dcrdp1x = {0, (p3 - p2).z(), -(p3 - p2).y()},
               dcrdp1y = {-(p3 - p2).z(), 0, (p3 - p2).x()},
               dcrdp1z = {(p3 - p2).y(), -(p3 - p2).x(), 0};
            TV dcrdp2x = {0, -(p3 - p1).z(), (p3 - p1).y()},
               dcrdp2y = {(p3 - p1).z(), 0, -(p3 - p1).x()},
               dcrdp2z = {-(p3 - p1).y(), (p3 - p1).x(), 0};
            TV dcrdp3x = {0, (p2 - p1).z(), -(p2 - p1).y()},
               dcrdp3y = {-(p2 - p1).z(), 0, (p2 - p1).x()},
               dcrdp3z = {(p2 - p1).y(), -(p2 - p1).x(), 0};

            if (cr.squaredNorm() == 0)
                continue;
            T normCoef = -pow(cr.squaredNorm(), -1.5);
            T normI = 1 / cr.norm();
            TV dnormdp1 = {p1.x() * ((p3 - p2).z() * (p3 - p2).z() +
                                     (p3 - p2).y() * (p3 - p2).y()),
                           p1.y() * ((p3 - p2).z() * (p3 - p2).z() +
                                     (p3 - p2).x() * (p3 - p2).x()),
                           p1.z() * ((p3 - p2).x() * (p3 - p2).x() +
                                     (p3 - p2).y() * (p3 - p2).y())};
            dnormdp1 *= normCoef;

            TV dnormdp2 = {p2.x() * ((p3 - p1).z() * (p3 - p1).z() +
                                     (p3 - p1).y() * (p3 - p1).y()),
                           p2.y() * ((p3 - p1).z() * (p3 - p1).z() +
                                     (p3 - p1).x() * (p3 - p1).x()),
                           p2.z() * ((p3 - p1).x() * (p3 - p1).x() +
                                     (p3 - p1).y() * (p3 - p1).y())};
            dnormdp2 *= normCoef;

            TV dnormdp3 = {p3.x() * ((p1 - p2).z() * (p1 - p2).z() +
                                     (p1 - p2).y() * (p1 - p2).y()),
                           p3.y() * ((p1 - p2).z() * (p1 - p2).z() +
                                     (p1 - p2).x() * (p1 - p2).x()),
                           p3.z() * ((p1 - p2).x() * (p1 - p2).x() +
                                     (p1 - p2).y() * (p1 - p2).y())};
            dnormdp3 *= normCoef;

            T d = fromBelow * dist - h;
            TV gradq, gradp1, gradp2, gradp3;
            gradq << fromBelow * n;
            gradp1 << fromBelow * (-n + dnormdp1 * ((q - p1).dot(cr)) +
                                   normI * TV((q - p1).dot(dcrdp1x),
                                              (q - p1).dot(dcrdp1y),
                                              (q - p1).dot(dcrdp1z)));
            gradp2 << fromBelow * (dnormdp2 * ((q - p1).dot(cr)) +
                                   normI * TV((q - p1).dot(dcrdp2x),
                                              (q - p1).dot(dcrdp2y),
                                              (q - p1).dot(dcrdp2z)));
            gradp3 << fromBelow * (dnormdp3 * ((q - p1).dot(cr)) +
                                   normI * TV((q - p1).dot(dcrdp3x),
                                              (q - p1).dot(dcrdp3y),
                                              (q - p1).dot(dcrdp3z)));

            T squaredGradientsSum =
                (gradq.squaredNorm() + gradp1.squaredNorm() +
                 gradp2.squaredNorm() + gradp3.squaredNorm());
            if (squaredGradientsSum == 0)
                continue;
            T s = d / squaredGradientsSum;
            // n in formula (8) is the number of points (involved in the
            // constraint, so 4 here), NOT the normal vector
            T dp = -4 * s / sum_w;

            //        std::cout<<"delta pos: "<<dp*w_q*constraint.gradq<<" "<<dp
            //        * w_p1 * constraint.gradp1<<" "<<dp * w_p2 *
            //        constraint.gradp2<<" "<<dp * w_p3 *
            //        constraint.gradp3<<std::endl;

            p.row(constraint.qIdx) += dp * w_q * gradq;
            p.row(t1) += dp * w_p1 * gradp1;
            p.row(t2) += dp * w_p2 * gradp2;
            p.row(t3) += dp * w_p3 * gradp3;

            constraint.n = n * fromBelow;

            //        std::cout<<constraint.d<<" "<<std::abs(dist)-h<<std::endl;
        }
    }
}



void PBD::applyFriction()
{
    for (auto& constraint : collisionConstraintsList)
    {
        if (constraint.n == TV(0, 0, 0))
            continue;
        
        TV vRelative = v.row(constraint.qIdx);   // Current velocity
        TV vNormal = vRelative.dot(constraint.n) * constraint.n; // Component along the normal
        TV vTangential = vRelative - vNormal; // Component perpendicular to the normal

        // Apply friction: scale tangential component
        TV vTangentialFriction = (1.0 - mu) * vTangential;

        // Update velocity
        v.row(constraint.qIdx) = vTangentialFriction + vNormal;
    }
}

void PBD::positionConstraints()
{
    for (int i = 0; i < posConstraintsIdxs.rows(); ++i)
    {
        int idx = posConstraintsIdxs[i];
        p.row(idx) = posConstraintsV.row(i);
    }
    if (floorCollision)
    {
        for (int i = 0; i < p.rows(); ++i)
        {
            // floor at z=-1/4
            p(i, 2) = std::max(p(i, 2), -0.25);
        }
    }
}

void PBD::projectConstraints(int solver_it)
{
    if (stretchingConstraintsActivated)
        if (useXPBD)
            PBD::stretchingConstraintsXPBD();
        else
            PBD::stretchingConstraints(solver_it);


    if (bendingConstraintsActivated)
        if (useXPBD)
            PBD::bendingConstraintsXPBD();
        else
            PBD::bendingConstraints(solver_it);

    if (collisionConstraintsActivated)
        PBD::collisionConstraints();

    if (positionConstraintsActivated)
        PBD::positionConstraints();
}

void PBD::dampVelocities(T kDamping)
{
    TV xCm = TV(0, 0, 0);
    TV vCm = TV(0, 0, 0);
    T mTot = m(0);
    for (int i = 0; i < currentV.rows(); i++)
    {
        xCm += currentV.row(i) * m(i);
        vCm += v.row(i) * m(i);
        mTot += m(i);
    }
    xCm /= mTot;
    vCm /= mTot;
    TV L = TV(0, 0, 0);
    TM I;
    I << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    MatrixXT r(currentV.rows(), 3);
    for (int i = 0; i < currentV.rows(); i++)
    {
        TV ri = TV(currentV.row(i)) - xCm;
        L += ri.cross(TV(v.row(i) * m(i)));
        TM ridash;
        ridash << 0, -ri.z(), ri.y(), ri.z(), 0, -ri.x(), -ri.y(), ri.x(), 0;
        I += ridash * ridash.transpose() * m(i);
        r.row(i) = ri;
    }
    TV omega = I.inverse() * L;
    for (int i = 0; i < currentV.rows(); i++)
    {
        TV deltaVi = vCm + omega.cross(TV(r.row(i))) - TV(v.row(i));
        v.row(i) += kDamping * deltaVi;
    }
}

void PBD::applyFakeWind()
{
    int num_faces = 3;
    for (int i = 0; i < num_faces; ++i)
    {
        int f = distF(gen);
        TV force = { -dist10(gen) * 2.0, dist10(gen) * 0.001 - 0.005, dist10(gen) * 0.01 - 0.05};
        for (int j = 0; j < 3; ++j)
        {
            v.row(faces(f, j)) += dt * force;
        }
    }
}


bool PBD::advanceOneStep(int step)
{
    //    std::cout<<step<<"\n";
    // Pre-solve: apply external forces
    int nRows = currentV.rows();
    for (int i = 0; i < nRows; i++)
    {
        v.row(i) += dt * g;

        if (fakeWindActivated)
        {
            int wind = dist10(gen);
            if (wind == 10) applyFakeWind();
        }
    }

    dampVelocities(k_damping);

    for (int i = 0; i < nRows; i++)
    {
        // Verlet integration, as in paper
        p.row(i) = currentV.row(i) + dt * v.row(i);
    }

    // for XPBD
    lambdas.setZero(constraint_idx.size());

    // generate collision constraints
    if (collisionConstraintsActivated)
        generateCollisionConstraints();

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
    applyFriction();

    if (step == nSteps)
        return true;
    return false;
}

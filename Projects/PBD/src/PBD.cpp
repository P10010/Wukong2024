#include "../include/PBD.h"
#include "igl/Hit.h"
#include <igl/ray_mesh_intersect.h>
#include <cmath>
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/tri_tri_intersect.h>
#include <random>
#include <iomanip>

void PBD::initializeFromFile(const std::string& filename)
{
    scene = filename.data();
    MatrixXT V;
    MatrixXi F;
    igl::readOBJ(scene, V, F);



    igl::readOBJ("../../../Projects/PBD/data/sailboat/parts/boat.obj", boatV, boatF);


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
//    std::cout<<nRows<<" "<<p.rows()<<std::endl;
    m.resize(nRows);
    w.resize(nRows);

    testPoint=currentV.row(0);

//  std::cout<<currentV.cols()<<std::endl;
//
//    for (int i = 0; i < nRows; i++) {
//      for (int j = 0; j < currentV.cols(); j++)
//        std::cout << currentV(i, j) << " ";
//      std::cout<<std::endl;
//    }



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

    adjFaces=std::vector<std::vector<int>>(nRows);

    for(int i=0;i<faces.rows();i++)
      for(int j=0;j<faces.cols();j++) {
        int vertex=faces(i,j);
        adjFaces[vertex].push_back(i);
      }

    igl::edges(F, edges);
    edge_lengths.resize(edges.rows());
    for (int i = 0; i < edges.rows(); i++)
    {
        int u = edges(i, 0);
        int v = edges(i, 1);
        edge_lengths(i) = (atRest.row(u) - atRest.row(v)).norm();
        maxEdgeLength=std::max(edge_lengths(i),maxEdgeLength);
    }

    igl::edges(boatF,boatEdges);
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

    currentBoatV=boatV;
    int nRowsBoat = boatV.rows();
    p.resize(nRowsBoat, 3);
    boatM.resize(nRowsBoat);
    boatW.resize(nRowsBoat);

    // Calculate mass for each vertex based on adjacent triangles
    for (int i = 0; i < nRowsBoat; i++)
    {
      boatM(i) = 10000;
    }
//    for (int i = 0; i < faces.rows(); i++)
//    {
//      int v1 = faces(i, 0);
//      int v2 = faces(i, 1);
//      int v3 = faces(i, 2);
//      TV a = atRest.row(v2) - atRest.row(v1);
//      TV b = atRest.row(v3) - atRest.row(v1);
//      T area = a.cross(b).norm() * 0.5;
//      m(v1) += area / 3.0;
//      m(v2) += area / 3.0;
//      m(v3) += area / 3.0;
//    }

    for (int i = 0; i < nRowsBoat; i++)
    {
      boatW(i) = 0 ;
    }

    v = MatrixXT::Zero(nRows, 3);
//    Pinned to diagonal and initialization for incorrect edge collisions
//    posConstraintsIdxs.resize(2);
//    posConstraintsIdxs << 9, 90;
//    posConstraintsV.resize(2, 3);
//    posConstraintsV << -0.3, 0.0, 0.0, 0.0, 0.3, 0.0;
//    w(9) = 0;
//    w(90) = 0;
//
//    incidentFaces.resize(edges.rows());
//    for (int e = 0; e < edges.rows(); e++)
//        for(int f=0; f< faces.rows();f++){
//            int cnt=0;
//            for(int i=0;i<2;i++)
//                for(int j=0;j<3;j++)
//                    cnt+=edges(e,i)==faces(f,j);
//            if(cnt==2)
//                incidentFaces[e].push_back(f);
//        }
//    adjList.reserve(atRest.rows());
//    for (int e = 0; e < edges.rows(); e++) {
//        adjList[edges(e, 0)].push_back(edges(e,1));
//        adjList[edges(e, 1)].push_back(edges(e,0));
//    }

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

//is q above the triangle (p1,p2,p3) with respect to the triangle normal?
int PBD::isAbove(const PBD::TV& q, const PBD::TV& p1, const PBD::TV& p2, const PBD::TV& p3) {
    // Compute triangle normal
    TV cr = (p2 - p1).cross(p3 - p1);
    TV n = cr.normalized();
    T dist = (q - p1).dot(n);

    return dist>0 ? 1 : -1;
}

//returns the squared norm
T PBD::closestToPointInTriangle(const PBD::TV& q, const PBD::TV& p1, const PBD::TV& p2, const PBD::TV& p3, PBD::TV& closestPoint){
  // Compute triangle normal
  TV cr = (p2 - p1).cross(p3 - p1);
  TV n = cr.normalized();
  T dist = (q - p1).dot(n);

  // calculate barycentric coordinates of projection into the plane (from
  // Robust Treatment of Collisions, Contact and Friction for Cloth Animation)
  TV x13 = p1 - p3, x23 = p2 - p3, x43 = q - p3;
  TM2 A;
  A << x13.dot(x13), x13.dot(x23), x13.dot(x23), x23.dot(x23);
  TV2 b = {x13.dot(x43), x23.dot(x43)};
  TV2 w12 = A.colPivHouseholderQr().solve(b);


//  //Clamp the barycentric coordinates to 0,1 so that they lie inside the triangle
//  std::clamp(w12[0],0.,1.);
//  std::clamp(w12[1],0.,1.);

  T w3 = 1 - w12[0] - w12[1];

  //closest point on the plane is not inside the triangle -> choose the vertex p1 (todo)
  if(w12[0]<0 || w12[0]>1 || w12[1]<0 || w12[1]>1 || w3<0 || w3>1)
    closestPoint=p1;
  else
    closestPoint=w12[0]*p1+w12[1]*p2+w3*p3;

  //return squared norm of the distance
  return (closestPoint-q).squaredNorm();
}

// check whether point q is closer than the given thickness to the triangle with vertices p1, p2, p3
bool PBD::pointIntersectsTriangle(const PBD::TV& q, const PBD::TV& p1,
                                  const PBD::TV& p2, const PBD::TV& p3, T thickness)
{
    // Compute triangle normal
    TV cr = (p2 - p1).cross(p3 - p1);
    TV n = cr.normalized();

    // Compute the distance of the point from the plane containing the triangle
    T dist = (q - p1).dot(n);

    // If the dist is greater than the thickness the point is not close to the plane
    // containing the triangle
    if (abs( dist) > thickness)
        return false;

    // calculate barycentric coordinates of projection into the plane (from
    // Robust Treatment of Collisions, Contact and Friction for Cloth Animation)
    TV x13 = p1 - p3, x23 = p2 - p3, x43 = q - p3;
    TM2 A;
    A << x13.dot(x13), x13.dot(x23), x13.dot(x23), x23.dot(x23);
    TV2 b = {x13.dot(x43), x23.dot(x43)};
    TV2 w12 = A.colPivHouseholderQr().solve(b);
    T w3 = 1 - w12[0] - w12[1];
    // δ is the thickness divided by a characteristic length of the triangle, i.e. squared
    // root of the area
    T area = 0.5 * cr.norm();
    if (area == 0)
        return false;
    T delta = thickness / (sqrt(abs(area)));
    return (w12.x() >= -delta && w12.x() <= 1 + delta && w12.y() >= -delta &&
            w12.y() <= 1 + delta && w3 >= -delta &&
            w3 <= 1 + delta); // collision detected
}

// check if edge x1x2 is close to x3x4 (from
// Robust Treatment of Collisions, Contact and Friction for Cloth Animation)
bool PBD::edgesAreClose(const PBD::TV& x1, const PBD::TV& x2, const PBD::TV& x3, const PBD::TV& x4) const{

  TV x21=x2-x1, x43=x4-x3;

  //parallel edges
  if(x21.cross(x43).norm()<epsilon) {
    //currently I only check the distance between the two lines, I should project the points to do a proper calculation

    //horizontal lines
    if(x21.x()<epsilon)
      return abs(x2.x()-x4.x())<h;

    float m=x21.y()/x21.x();
    float q1=x2.y()-m*x2.x();
    float q2=x4.y()-m*x4.x();

    return abs(q1-q2)/sqrt(m*m+1)<h;
  }

  //Find the closest points on the infinite lines
  TV x31=x3-x1;
  TM2 A;
  A << x21.dot(x21), -x21.dot(x43), -x21.dot(x43), x43.dot(x43);
  TV2 b = {x21.dot(x31), -x43.dot(x31)};

  TV2 w12 = A.colPivHouseholderQr().solve(b);

  if(w12.x()>1)
    w12.x()=1;
  if(w12.y()>1)
    w12.y()=1;
  if(w12.x()<0)
    w12.x()=0;
  if(w12.y()<0)
    w12.y()=0;
  TV point1=x1+w12.x()*x21, point2=x3+w12.y()*x43;

  return (point1-point2).norm()<h;
}

int PBD::hash(int i, int j, int k, int n)
{
    long hashedValue=(((long)i * prime1) ^ ((long)j * prime2) ^ ((long)k * prime3));
    hashedValue%=n;
    return hashedValue;
}
int PBD::hash(const PBD::TV point, T l, int n, const TV& minCoord)
{
    int i = std::floor((point.x() - minCoord.x()) / l),
        j = std::floor((point.y() - minCoord.y()) / l),
        k = std::floor((point.z() - minCoord.z()) / l);
//    std::cout << point.x() << " " << point.y() << point.z() << " " << i << " " << j << " " << k << " " << n << " " << hash(i, j, k, n) << std::endl;
    return hash(i, j, k, n);
}

void PBD::hashVertices(std::vector<std::vector<int>>& hashTable, T boxSize,
                       TV& minCoord)
{
    int n = hashTable.size();
//    if(p.rows()!=currentV.rows())
//    std::cout<<p.rows()<<" "<<currentV.rows()<<"\n";
    int nRows=currentV.rows();
    for (int i = 0; i < nRows; i++)
    {
        TV q = p.row(i);
//        std::cout<<std::setprecision(20);
//        std::cout<<i<<" "<<q.x()<<" "<<q.y()<<" "<<q.z()<<std::endl;
        int hashVal = hash(q, boxSize, n, minCoord);
        hashTable[hashVal].push_back(i);
//        if(hashTable[hashVal].size()>10)
//          std::cout<<hashVal<<" "<<hashTable[hashVal].size()<<"\n";
    }

//    int boatRows=boatV.rows();
//    for (int i = 0; i < boatRows; i++)
//    {
//        TV q = boatV.row(i);
//        int hashVal = hash(q, boxSize, n, minCoord);
//        hashTable[hashVal].push_back(i+nRows);
//    }
}

void PBD::faceSelfCollisionConstraint(std::vector<std::vector<int>>& hashTable, const TV& minCoord, const TV& maxCoord,  int n, const T boxSize){
  for (int f = 0; f < faces.rows(); f++) {
    int t1 = faces(f, 0), t2 = faces(f, 1), t3 = faces(f, 2);
    TM triangleVertices;
    TV p1 = p.row(t1), p2 = p.row(t2), p3 = p.row(t3);
    triangleVertices << p1, p2, p3;
    std::array<int, 3> minBox, maxBox;
    for (int i = 0; i < 3; i++) {
      TV tmp = triangleVertices.row(i);
      minBox[i] = (int) std::floor(
              (triangleVertices.row(i).minCoeff() - minCoord(i)) / boxSize);
      maxBox[i] = (int) std::floor(
              (triangleVertices.row(i).maxCoeff() - minCoord(i)) / boxSize);
    }

    for (int i = minBox[0]; i <= maxBox[0]; i++)
      for (int j = minBox[1]; j <= maxBox[1]; j++)
        for (int k = minBox[2]; k <= maxBox[2]; k++) {
          int hashVal = hash(i, j, k, n);
          for (int qIdx: hashTable[hashVal])
            if(qIdx!=t1 && qIdx!=t2 && qIdx!=t3)
              for(int f1 : adjFaces[qIdx]){
                int idx1=faces(f1,0), idx2=faces(f1,1), idx3=faces(f1,2);
                bool adj=false;
                for(int j=0;j<3;j++)
                  if(TT(f,j)==f1)
                    adj=true;
                if(!adj) {
                  Eigen::Matrix<double, 1, 3> r1 = p.row(idx1).transpose(), r2 = p.row(
                          faces(idx2)).transpose(), r3 = p.row(idx3).transpose();
                  bool coplanar;
                  Eigen::Matrix<double, 1, 3> src, target;
                  Eigen::Matrix<double, 1, 3> p1T = p1.transpose(), p2T = p2.transpose(), p3T = p3.transpose();
                  if (igl::tri_tri_intersection_test_3d(p1T, p2T, p3T, r1, r2, r3, coplanar, src, target))
                  if(!coplanar){
                    std::cout<<f<<" "<<f1<<std::endl;
                    CollisionConstraint constraint;
                    int qAbove = isAbove(p.row(qIdx), p1, p2, p3);
                    constraint.qIdx = qIdx;
                    constraint.f = faces.row(f);
                    constraint.inBoat = false;
                    constraint.fromBelow = -qAbove;
                    constraint.n = ((p2 - p1).cross(p3 - p1)).normalized();
                    collisionConstraintsList.push_back(constraint);
                    break;
                  }
                }
              }
        }
  }
}
void PBD::spatialHashing()
{
    int nRows=currentV.rows();
    TV minCoord = p.colwise().minCoeff();
    TV maxCoord = p.colwise().maxCoeff();
//    TV minCoordBoat = boatV.colwise().minCoeff();
//    TV maxCoordBoat = boatV.colwise().maxCoeff();
//    for(int i=0;i<3;i++) {
//        minCoord.coeffRef(i) = std::min(minCoord(i), minCoordBoat(i));
//        maxCoord.coeffRef(i) = std::max(maxCoord(i), maxCoordBoat(i));
//    }

    // T coeff = 1; //(pow(p.rows(),T(0.33)));
    T boxSize = ((maxCoord - minCoord).maxCoeff()) / 100;
//     T boxSize = 2 * h;
    int n = 4 * (p.rows());
    std::vector<std::vector<int>> hashTable(n);
    hashVertices(hashTable, boxSize, minCoord);
//    std::cout<<p.rows()<<std::endl;
//    faceSelfCollisionConstraint(hashTable,minCoord,maxCoord,n,boxSize);
//    return;
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
                    for (int qIdx : hashTable[hashVal]) {
                        //necessary condition for intersection: x_i -> p_i passes through the face
                        //This is equivalent to p_i being above the triangle and x_i below, or viceversa
                        //x_i is stored in currentV
                        if(qIdx>nRows){
                            qIdx-=nRows;
                            //todo fix this
//                            if (pointIntersectsTriangle(boatV.row(qIdx), p1, p2, p3)) {
//                                CollisionConstraint constraint;
//                                constraint.qIdx = qIdx;
//                                constraint.f = faces.row(f);
//                                constraint.fromBelow = -isAbove(boatV.row(qIdx),p1,p2,p3);
//                                constraint.inBoat = false;
//                                collisionConstraintsList.push_back(constraint);
//                            }
                        }
                        else {
                            int qAbove = isAbove(p.row(qIdx), p1, p2, p3);
                            if (qAbove != isAbove(currentV.row(qIdx), p1, p2, p3) &&
                                qIdx != t1 && qIdx != t2 && qIdx != t3 &&
                                pointIntersectsTriangle(p.row(qIdx), p1, p2, p3,h)) {
                                CollisionConstraint constraint;
                                constraint.qIdx = qIdx;
                                constraint.f = faces.row(f);
                                constraint.fromBelow = -qAbove;
                                constraint.inBoat = false;
                                collisionConstraintsList.push_back(constraint);
                            }
                        }
                    }
                }
    }


//    for (int f = 0; f < boatF.rows(); f++)
//    {
//      int t1 = boatF(f, 0), t2 = boatF(f, 1), t3 = boatF(f, 2);
//      TM triangleVertices;
//      TV p1 = boatV.row(t1), p2 = boatV.row(t2), p3 = boatV.row(t3);
//      triangleVertices << p1, p2, p3;
//      std::array<int, 3> minBox, maxBox;
//      // std::cout << "p1\n" << p1 << std::endl;
//      // std::cout << "triangle vertices\n"
//      //           << triangleVertices << std::endl;
//      for (int i = 0; i < 3; i++)
//      {
//        TV tmp = triangleVertices.row(i);
//        minBox[i] = (int)std::floor(
//                (triangleVertices.row(i).minCoeff() - minCoord(i)) / boxSize);
//        if(minBox[i]<0)
//            minBox[i]=0;
//
//        if(triangleVertices.row(i).maxCoeff()>maxCoord(i))
//            maxBox[i] = (int)std::floor(
//                    (maxCoord(i) - minCoord(i)) / boxSize);
//        else
//            maxBox[i] = (int)std::floor(
//                    (triangleVertices.row(i).maxCoeff() - minCoord(i)) / boxSize);
//        //todo fix this properly
//        if(maxBox[i]<minBox[i])
//            maxBox[i]=minBox[i];
//      }
////       std::cout << "minBox\n" << minBox[0] << " " << minBox[1] << " "
////                 << minBox[2] << std::endl;
////       std::cout << "maxBox\n" << maxBox[0] << " " << maxBox[1] << " "
////                   << maxBox[2] << std::endl;
//
//      for (int i = minBox[0]; i <= maxBox[0]; i++)
//        for (int j = minBox[1]; j <= maxBox[1]; j++)
//          for (int k = minBox[2]; k <= maxBox[2]; k++)
//          {
//            int hashVal = hash(i, j, k, n);
//            for (int qIdx : hashTable[hashVal]) {
//                //necessary condition for intersection: x_i -> p_i passes through the face
//                //This is equivalent to p_i being above the triangle and x_i below, or viceversa
//                int qAbove=isAbove(p.row(qIdx), p1, p2, p3);
////                if(p.row(qIdx).z()<triangleVertices.row(2).minCoeff())
//                if (qAbove != isAbove(currentV.row(qIdx), p1, p2, p3)
//                    && pointIntersectsTriangle(p.row(qIdx), p1, p2, p3,maxEdgeLength))
//                {
////                    if (w(qIdx) > 0)
////                        std::cout << "its with boat detected" << qIdx << " " << f << std::endl;
//                    CollisionConstraint constraint;
//                    constraint.qIdx = qIdx;
//                    constraint.f = boatF.row(f);
//                    constraint.inBoat = true;
//                    constraint.fromBelow = -qAbove;
//                    collisionConstraintsList.push_back(constraint);
//                }
//            }
//          }
//    }
// wrong implementation of edge-edge collisions
//  for (int e = 0; e < edges.rows(); e++)
//  {
//    int a, b;
//    a = edges(e, 0);
//    b = edges(e, 1);
//    VMD::Matrix<T, 3, 2> edgeVertices;
//    TV p1=p.row(a), p2=p.row(b);
//    edgeVertices << p1,p2;
//    std::array<int, 3> minBox, maxBox;
//    // std::cout << "p1\n" << p1 << std::endl;
//    // std::cout << "triangle vertices\n"
//    //           << triangleVertices << std::endl;
//    for (int i = 0; i < 3; i++)
//    {
//      TV2 tmp = edgeVertices.row(i);
//      minBox[i] = (int)std::floor(
//              (tmp.minCoeff() - minCoord(i)) / boxSize);
//      // std::max(0,(int)std::floor((triangleVertices.row(i).minCoeff()-minCoord(i))/boxSize)-1);
//      maxBox[i] = (int)std::floor(
//              (tmp.maxCoeff() - minCoord(i)) / boxSize);
//    }
////     std::cout << "minBox\n" << minBox[0] << " " << minBox[1] << " "
////               << minBox[2] << std::endl;
////     std::cout << "maxBox\n" << maxBox[0] << " " << maxBox[1] << " "
////                 << maxBox[2] << std::endl;
//
//    for (int i = minBox[0]; i <= maxBox[0]; i++)
//      for (int j = minBox[1]; j <= maxBox[1]; j++)
//        for (int k = minBox[2]; k <= maxBox[2]; k++)
//        {
//          int hashVal = hash(i, j, k, n);
//          for (int qIdx : hashTable[hashVal])
//            for(int endPoint : adjList[qIdx])
//              if (qIdx != a && qIdx != b && endPoint!=a && endPoint!=b &&
//                  edgesAreClose(p.row(qIdx), p.row(endPoint), p1, p2))
//              {
//                CollisionConstraint constraint;
//                constraint.qIdx = qIdx;
//                constraint.f = faces.row(incidentFaces[e][0]);
//                collisionConstraintsList.push_back(constraint);
//
////                constraint.qIdx = endPoint;
////                collisionConstraintsList.push_back(constraint);
//              }
//        }
//  }
}

void PBD::faceStaticConstraint(std::vector<std::vector<int>>& hashTable, const TV& minCoord, const TV& maxCoord,  int n, const T boxSize){
  for (int f = 0; f < boatF.rows(); f++) {
    int t1 = boatF(f, 0), t2 = boatF(f, 1), t3 = boatF(f, 2);
    TM triangleVertices;
    TV p1 = boatV.row(t1), p2 = boatV.row(t2), p3 = boatV.row(t3);
    triangleVertices << p1, p2, p3;
    std::array<int, 3> minBox{}, maxBox{};
    for (int i = 0; i < 3; i++) {
      TV tmp = triangleVertices.row(i);
      minBox[i] = (int) std::floor(
              (triangleVertices.row(i).minCoeff() - minCoord(i)) / boxSize);
      if (minBox[i] < 0)
        minBox[i] = 0;

      if (triangleVertices.row(i).maxCoeff() > maxCoord(i))
        maxBox[i] = (int) std::floor(
                (maxCoord(i) - minCoord(i)) / boxSize);
      else
        maxBox[i] = (int) std::floor(
                (triangleVertices.row(i).maxCoeff() - minCoord(i)) / boxSize);
      //todo fix this properly
      if (maxBox[i] < minBox[i])
        maxBox[i] = minBox[i];
    }

    for (int i = minBox[0]; i <= maxBox[0]; i++)
      for (int j = minBox[1]; j <= maxBox[1]; j++)
        for (int k = minBox[2]; k <= maxBox[2]; k++) {
          int hashVal = hash(i, j, k, n);
          for (int qIdx: hashTable[hashVal])
            for(int f1 : adjFaces[qIdx]){
              Eigen::Matrix<double, 1, 3> r1 = p.row(faces(f1,0)).transpose(), r2 = p.row(faces(f1,1)).transpose(), r3 = p.row(faces(f1,2)).transpose();
              bool coplanar;
              Eigen::Matrix<double, 1, 3> src,target;
              Eigen::Matrix<double, 1, 3> p1T=p1.transpose(), p2T=p2.transpose(), p3T=p3.transpose();
              if(igl::tri_tri_intersection_test_3d(p1T,p2T,p3T,r1,r2,r3,coplanar,src,target)){
                CollisionConstraintStatic constraint;
                constraint.pIdx = qIdx;
                constraint.q= src;
                constraint.n = ((p2 - p1).cross(p3 - p1)).normalized();
                collisionConstraintsStaticList.push_back(constraint);
                break;
              }
            }
        }
  }
}
void PBD::edgeStaticConstraint(std::vector<std::vector<int>>& hashTable, const TV& minCoord, const TV& maxCoord, int n, const T boxSize){
  for (int e = 0; e < boatEdges.rows(); e++)
  {
    int a, b;
    a = edges(e, 0);
    b = edges(e, 1);
    VMD::Matrix<T, 3, 2> edgeVertices;
    TV p1=p.row(a), p2=p.row(b);
    edgeVertices << p1,p2;
    std::array<int, 3> minBox{}, maxBox{};

    for (int i = 0; i < 3; i++)
    {
      TV2 tmp = edgeVertices.row(i);
      minBox[i] = (int)std::floor(
              (edgeVertices.row(i).minCoeff() - minCoord(i)) / boxSize);
      if(minBox[i]<0)
        minBox[i]=0;

      if(edgeVertices.row(i).maxCoeff()>maxCoord(i))
        maxBox[i] = (int)std::floor(
                (maxCoord(i) - minCoord(i)) / boxSize);
      else
        maxBox[i] = (int)std::floor(
                (edgeVertices.row(i).maxCoeff() - minCoord(i)) / boxSize);
      //todo fix this properly
      if(maxBox[i]<minBox[i])
        maxBox[i]=minBox[i];
    }
//     std::cout << "minBox\n" << minBox[0] << " " << minBox[1] << " "
//               << minBox[2] << std::endl;
//     std::cout << "maxBox\n" << maxBox[0] << " " << maxBox[1] << " "
//                 << maxBox[2] << std::endl;

    for (int i = minBox[0]; i <= maxBox[0]; i++)
      for (int j = minBox[1]; j <= maxBox[1]; j++)
        for (int k = minBox[2]; k <= maxBox[2]; k++)
        {
          int hashVal = hash(i, j, k, n);
          for (int qIdx : hashTable[hashVal])
            for(int endPoint : adjList[qIdx])
              if (qIdx != a && qIdx != b && endPoint!=a && endPoint!=b &&
                  edgesAreClose(p.row(qIdx), p.row(endPoint), p1, p2))
              {
                CollisionConstraintStatic constraint;
                constraint.pIdx = qIdx;
                int f=incidentFaces[e][0];
                int t1 = boatF(f, 0), t2 = boatF(f, 1), t3 = boatF(f, 2);
                TM triangleVertices;
                TV p1 = boatV.row(t1), p2 = boatV.row(t2), p3 = boatV.row(t3);
//                constraint.n=
                collisionConstraintsStaticList.push_back(constraint);

                constraint.pIdx = endPoint;
                collisionConstraintsStaticList.push_back(constraint);
              }
        }
  }
}
void PBD::spatialHashingStatic() {
  int nRows = currentV.rows();
  TV minCoord = p.colwise().minCoeff();
  TV maxCoord = p.colwise().maxCoeff();
//    TV minCoordBoat = boatV.colwise().minCoeff();
//    TV maxCoordBoat = boatV.colwise().maxCoeff();
//    for(int i=0;i<3;i++) {
//        minCoord.coeffRef(i) = std::min(minCoord(i), minCoordBoat(i));
//        maxCoord.coeffRef(i) = std::max(maxCoord(i), maxCoordBoat(i));
//    }

  // T coeff = 1; //(pow(p.rows(),T(0.33)));
//  T boxSize = ((maxCoord - minCoord).maxCoeff()) / 100;
  T boxSize = maxEdgeLength;
  int n = 4 * (p.rows());
  std::vector<std::vector<int>> hashTable(n);
  hashVertices(hashTable, boxSize, minCoord);

  faceStaticConstraint(hashTable,minCoord,maxCoord,n,boxSize);
  return;

  //for all the segments xi->pi that lie completely inside the mesh, find the closest point on the mesh
  std::vector<CollisionConstraintStatic> closestPoints(p.rows());
  //-1: not checked; 0: pi is inside, but xi is outside; 1: segment completely inside the mesh
  std::vector<int> rayInsideMesh(p.rows(),-1);

  std::vector<bool> toCheck(p.rows(),false);
  for (int f = 0; f < boatF.rows(); f++)
    {
      int t1 = boatF(f, 0), t2 = boatF(f, 1), t3 = boatF(f, 2);
      TM triangleVertices;
      TV p1 = boatV.row(t1), p2 = boatV.row(t2), p3 = boatV.row(t3);
      triangleVertices << p1, p2, p3;
      std::array<int, 3> minBox{}, maxBox{};
      // std::cout << "p1\n" << p1 << std::endl;
      // std::cout << "triangle vertices\n"
      //           << triangleVertices << std::endl;
      for (int i = 0; i < 3; i++)
      {
        TV tmp = triangleVertices.row(i);
        minBox[i] = (int)std::floor(
                (triangleVertices.row(i).minCoeff() - minCoord(i)) / boxSize);
        if(minBox[i]<0)
            minBox[i]=0;

        if(triangleVertices.row(i).maxCoeff()>maxCoord(i))
            maxBox[i] = (int)std::floor(
                    (maxCoord(i) - minCoord(i)) / boxSize);
        else
            maxBox[i] = (int)std::floor(
                    (triangleVertices.row(i).maxCoeff() - minCoord(i)) / boxSize);
        //todo fix this properly
        if(maxBox[i]<minBox[i])
            maxBox[i]=minBox[i];
      }
//       std::cout << "minBox\n" << minBox[0] << " " << minBox[1] << " "
//                 << minBox[2] << std::endl;
//       std::cout << "maxBox\n" << maxBox[0] << " " << maxBox[1] << " "
//                   << maxBox[2] << std::endl;

      for (int i = minBox[0]; i <= maxBox[0]; i++)
        for (int j = minBox[1]; j <= maxBox[1]; j++)
          for (int k = minBox[2]; k <= maxBox[2]; k++)
          {
            int hashVal = hash(i, j, k, n);
            for (int qIdx : hashTable[hashVal]) {
                TV closestPoint;
                TV xi=currentV.row(qIdx);
                TV pi=p.row(qIdx);
                T distSqr=closestToPointInTriangle(pi,p1,p2,p3,closestPoint);
                if(!toCheck[qIdx] || distSqr<(pi-closestPoints[qIdx].q).squaredNorm()) {
  //                        std::cout<<qIdx<<"  "<<closestPoint.x()<<" "<<closestPoint.y()<<" "<<closestPoint.z()<<" "<<distSqr<<std::endl;
                  closestPoints[qIdx].pIdx = qIdx;
                  closestPoints[qIdx].q=closestPoint;
                  closestPoints[qIdx].n = ((p2 - p1).cross(p3 - p1)).normalized();
                }
                toCheck[qIdx]=true;

                continue;
                //necessary condition for intersection: x_i -> p_i passes through the face
                //A necessary condition far that is p_i being above the triangle and x_i below, or viceversa
                int qAbove=isAbove(p.row(qIdx), p1, p2, p3);

                TV rayDir=pi-xi;
                igl::Hit hit{};
//                if(p.row(qIdx).z()<triangleVertices.row(2).minCoeff())
//qAbove != isAbove(currentV.row(qIdx), p1, p2, p3) && pointIntersectsTriangle(p.row(qIdx), p1, p2, p3,maxEdgeLength)
                if (igl::ray_triangle_intersect(xi,rayDir,boatV,boatF,f,hit))
//                  if(hit.t<1.f)
                {
//                    if (w(qIdx) > 0)
//                        std::cout << "its with boat detected" << qIdx << " " << f << std::endl;
                    bool firstIts;
                    if(rayInsideMesh[qIdx]==-1){
                      std::vector<igl::Hit> hits;
                      //this is slow, but hopefully this check will not be performed too often...
                      //it can probably be replaced by a clever ray-triangle intersection test... or maybe just doing spatial hashing for xi too
                      igl::ray_mesh_intersect(xi,rayDir,boatV,boatF,hits);
                      rayInsideMesh[qIdx]=hits.size()%2; //==1 ? 0 : 1;
                      firstIts=true;
//                      std::cout<<qIdx<<" "<<pi<<" "<<hits.size()<<"\n";
                    }
                    else
                      firstIts=false;
                    if( rayInsideMesh[qIdx]==0 ) {
                      if(hit.t<1.f) {
                        CollisionConstraintStatic constraint;
                        constraint.pIdx = qIdx;
                        constraint.q = xi + rayDir * hit.t;
                        constraint.n = ((p2 - p1).cross(p3 - p1)).normalized();
                        collisionConstraintsStaticList.push_back(constraint);
                      }
                    }
                    else{
                      TV closestPoint;
                      T distSqr=closestToPointInTriangle(pi,p1,p2,p3,closestPoint);
                      if(firstIts || distSqr<(pi-closestPoints[qIdx].q).squaredNorm()) {
//                        std::cout<<qIdx<<"  "<<closestPoint.x()<<" "<<closestPoint.y()<<" "<<closestPoint.z()<<" "<<distSqr<<std::endl;
                        closestPoints[qIdx].pIdx = qIdx;
                        closestPoints[qIdx].q=closestPoint;
                        closestPoints[qIdx].n = ((p2 - p1).cross(p3 - p1)).normalized();
                      }
                    }
                }
            }
          }
    }

  for(int i=0;i<p.rows();i++)
    if(rayInsideMesh[i]==1) {
      collisionConstraintsStaticList.push_back(closestPoints[i]);
//      std::cout<<i<<" ";
//      for(int j=0;j<3;j++)
//       std::cout<<closestPoints[i].q[j]<<" ";
//      std::cout<<std::endl;
    }

  for (int i = 0; i < p.rows(); ++i) {
    if(toCheck[i]){
      //check if the ray x_i -> p_i intersects or is inside the boat
      TV rayDir = p.row(i) - currentV.row(i);
      std::vector<igl::Hit> hits;
      TV xi = currentV.row(i);
      TV pi=p.row(i);
      if (igl::ray_mesh_intersect(xi, rayDir, boatV, boatF, hits)) {
        bool found = false;
        if (hits.size() % 2 == 0) {
          for (auto hit: hits) {
            if (hit.t < 1.f) {
  //          std::cout<<i<<" "<<hit.t<<std::endl;
              CollisionConstraintStatic constraint;
              constraint.pIdx = i;
              TV tmp = hit.t * rayDir;
              constraint.q = xi + tmp;
              TV p1 = boatV.row(boatF(hit.id, 0)), p2 = boatV.row(boatF(hit.id, 1)), p3 = boatV.row(boatF(hit.id, 2));
              constraint.n = ((p2 - p1).cross(p3 - p1)).normalized();
              collisionConstraintsStaticList.push_back(constraint);
              found = true;
              break;
            }
          }
        } else {


          TV zAxis = {0,0,1};
          if (!found && igl::ray_mesh_intersect(pi, zAxis, boatV, boatF, hits))
            if (hits[0].t <= maxEdgeLength) {
  //            std::cout << i << std::endl;
              auto hit = hits[0];
              CollisionConstraintStatic constraint;
              constraint.pIdx = i;
              TV tmp = hit.t * zAxis;
              constraint.q = pi + tmp;
              constraint.n = zAxis;
              collisionConstraintsStaticList.push_back(constraint);
            }

        }
        /*
        if(!found && hits.size()%2!=0)
          std::cout<<i<<" "<<hits.size()<<std::endl;
        if (!found && igl::ray_mesh_intersect(xi, -rayDir, boatV, boatF, hits)) {
          std::cout<<i<<std::endl;
          auto hit = hits[0];
          CollisionConstraintStatic constraint;
          constraint.pIdx = i;
          TV tmp = -hit.t * rayDir;
          constraint.q = xi + tmp;
          TV p1 = boatV.row(boatF(hit.id, 0)), p2 = boatV.row(boatF(hit.id, 1)), p3 = boatV.row(boatF(hit.id, 2));
          constraint.n = ((p2 - p1).cross(p3 - p1)).normalized();
          collisionConstraintsStaticList.push_back(constraint);
        }
         */
      }
    }
  }

//  std::cout<<collisionConstraintsStaticList.size()<<std::endl;
}
void PBD::generateCollisionConstraintsStatic(){
  // Clear previous collision constraints with static objects
  collisionConstraintsStaticList.clear();

  if(useSpatialHashing)
    spatialHashingStatic();
  else {
    // Iterate through all vertices
    for (int i = 0; i < p.rows(); ++i) {
      //check if the ray x_i -> p_i intersects or is inside the boat
      TV rayDir = p.row(i) - currentV.row(i);
      std::vector<igl::Hit> hits;
      TV xi = currentV.row(i);
      TV pi=p.row(i);
      if (igl::ray_mesh_intersect(xi, rayDir, boatV, boatF, hits)) {
        bool found = false;
        if(hits.size()%2==0) {
          for (auto hit: hits) {
            if (hit.t < 1.f) {
//          std::cout<<i<<" "<<hit.t<<std::endl;
              CollisionConstraintStatic constraint;
              constraint.pIdx = i;
              TV tmp = hit.t * rayDir;
              constraint.q = xi + tmp;
              TV p1 = boatV.row(boatF(hit.id, 0)), p2 = boatV.row(boatF(hit.id, 1)), p3 = boatV.row(boatF(hit.id, 2));
              constraint.n = ((p2 - p1).cross(p3 - p1)).normalized();
              collisionConstraintsStaticList.push_back(constraint);
              found = true;
              break;
            }
          }
        }
        else {
          TV zAxis = {0, 0, 1};
          if (!found && igl::ray_mesh_intersect(pi, zAxis, boatV, boatF, hits))
          if(hits[0].t<=maxEdgeLength){
//            std::cout << i << std::endl;
            auto hit = hits[0];
            CollisionConstraintStatic constraint;
            constraint.pIdx = i;
            TV tmp = hit.t * zAxis;
            constraint.q = pi + tmp;
            constraint.n = zAxis;
            collisionConstraintsStaticList.push_back(constraint);
          }
        }
        /*
        if(!found && hits.size()%2!=0)
          std::cout<<i<<" "<<hits.size()<<std::endl;
        if (!found && igl::ray_mesh_intersect(xi, -rayDir, boatV, boatF, hits)) {
          std::cout<<i<<std::endl;
          auto hit = hits[0];
          CollisionConstraintStatic constraint;
          constraint.pIdx = i;
          TV tmp = -hit.t * rayDir;
          constraint.q = xi + tmp;
          TV p1 = boatV.row(boatF(hit.id, 0)), p2 = boatV.row(boatF(hit.id, 1)), p3 = boatV.row(boatF(hit.id, 2));
          constraint.n = ((p2 - p1).cross(p3 - p1)).normalized();
          collisionConstraintsStaticList.push_back(constraint);
        }
         */
      }
    }

  }
  std::cout<<p.rows()<<" "<<collisionConstraintsStaticList.size()<<std::endl;

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
                int qAbove= isAbove(q, p1, p2, p3);
                if (qAbove != isAbove(currentV.row(i), p1, p2, p3) && pointIntersectsTriangle(q, p1, p2, p3,h))
                {
                    CollisionConstraint constraint;
                    constraint.qIdx = i;
                    constraint.f = faces.row(f);
                    constraint.fromBelow=-qAbove;
                    constraint.inBoat=false;
                    collisionConstraintsList.push_back(constraint);
                }
            }

//            for (int f = 0; f < boatF.rows(); f++)
//            {
//                int t1 = boatF(f, 0), t2 = boatF(f, 1), t3 = boatF(f, 2);
//                TV q = p.row(i);
//                TV p1 = boatV.row(t1), p2 = boatV.row(t2), p3 = boatV.row(t3);
//                int qAbove= isAbove(q, p1, p2, p3);
//                if (qAbove != isAbove(currentV.row(i), p1, p2, p3) && pointIntersectsTriangle(q, p1, p2, p3,maxEdgeLength))
//                {
//                    CollisionConstraint constraint;
//                    constraint.qIdx = i;
//                    constraint.f = boatF.row(f);
//                    constraint.fromBelow=-qAbove;
//                    constraint.inBoat=true;
//                    collisionConstraintsList.push_back(constraint);
//                }
//            }
        }
    }
//    std::cout<<"Number of dynamic collision constraints: "<<collisionConstraintsList.size()<<std::endl;
}

void PBD::collisionConstraintsStatic(){
  for (auto& constraint : collisionConstraintsStaticList)
  {
    TV itsPoint=p.row(constraint.pIdx);
    T dist=(itsPoint-constraint.q).dot(constraint.n);
//    std::cout<<dist<<std::endl;
    if(dist<0){
//      std::cout<<"idx: "<<constraint.pIdx<<" dist:"<<dist<<std::endl;
      TV grad=constraint.n;
      //only a point is involved and the squared norm of the normalized normal is 1
      p.row(constraint.pIdx)-=(dist)*grad;
    }
  }
}
void PBD::collisionConstraints()
{
    for (auto& constraint : collisionConstraintsList)
    {
        int t1 = constraint.f(0), t2 = constraint.f(1), t3 = constraint.f(2);

        TV q = p.row(constraint.qIdx);
        TV p1,p2,p3;
        if(!constraint.inBoat)
         p1 = p.row(t1), p2 = p.row(t2), p3 = p.row(t3);
        else
          p1 = boatV.row(t1), p2=boatV.row(t2), p3=boatV.row(t3);

        // Project the vertex and triangle vertices to satisfy the constraint
        T w_q = w(constraint.qIdx);
        T w_p1,w_p2,w_p3;
        if(!constraint.inBoat)
          w_p1 = w(t1), w_p2 = w(t2), w_p3 = w(t3);
        else
          w_p1 = boatW(t1), w_p2 = boatW(t2), w_p3 = boatW(t3);

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
//        if (dist < 0)
//            fromBelow = -1;
//        else
//            fromBelow = 1;
        fromBelow=constraint.fromBelow;
        if (fromBelow * dist < (constraint.inBoat ? maxEdgeLength : h))
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
            if(!constraint.inBoat) {
              p.row(t1) += dp * w_p1 * gradp1;
              p.row(t2) += dp * w_p2 * gradp2;
              p.row(t3) += dp * w_p3 * gradp3;
            }
            else{
              p.row(constraint.qIdx) += 10* dp * w_q * gradq;
//              std::cout<<"Collision with boat applied"<<constraint.qIdx<<" "<<constraint.f<<std::endl;
            }

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

    if (collisionConstraintsActivated){
      PBD::collisionConstraints();
      PBD::collisionConstraintsStatic();
    }


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
//    constantVelocity={constantXVelocity,0,0};
    // Pre-solve: apply external forces
    int nRows = currentV.rows();
    for (int i = 0; i < nRows; i++)
    {
        v.row(i) += dt * g;
        v.row(i).x() += constantXVelocity;

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

    for(int i=0;i<boatV.rows();i++)
      boatV.row(i).x()+=dt*constantXVelocity;

    // for XPBD
    lambdas.setZero(constraint_idx.size());

    // generate collision constraints
    if (collisionConstraintsActivated) {
      generateCollisionConstraints();
      generateCollisionConstraintsStatic();
    }

    for (int i = 0; i < posConstraintsIdxs.rows(); ++i) {
      int idx = posConstraintsIdxs[i];
      T dp = constantXVelocity * dt;
      posConstraintsV.row(i).x() += dp;
    }
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

    //velocity update
    applyFriction();

    if (step == nSteps)
        return true;
    return false;
}

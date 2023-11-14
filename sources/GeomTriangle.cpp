/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"

GeomTriangle::GeomTriangle() {
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy) {
    fNodeIndices = copy.fNodeIndices;

}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    fNodeIndices = copy.fNodeIndices;

    return *this;
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();

    phi[0] = 1.0 - xi[0] - xi[1];
    dphi(0, 0) = -1.0;
    dphi(1, 0) = -1.0;

    phi[1] = xi[0];
    dphi(0, 1) = 1.0;
    dphi(1, 1) = 0.0;

    phi[2] = xi[1];
    dphi(0, 2) = 0.0;
    dphi(1, 2) = 1.0;
}

void GeomTriangle::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);

    Shape(xi, phi, dphi);
    int space = NodeCo.rows();
    x.setZero();

    for (int i = 0; i < space; i++) {
        for (int j = 0; j < nCorners; j++){
            x[i] += NodeCo(i,j)*phi[j];
        }
    }
}

void GeomTriangle::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();  /// Duvida: por que NodeCo pode ter mais linhas que Dimensao do problema?
    if(NodeCo.cols() != nCorners) DebugStop();

    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);

    Shape(xi, phi, dphi);
    int space = NodeCo.rows();

    x.setZero();
    gradx.resize(space, Dimension);
    gradx.setZero();

    for (int i = 0; i < space; i++) {
        for (int j = 0; j < nCorners; j++) {
            x[i] += NodeCo(i,j) * phi[j];
            for (int k = 0; k < Dimension; k++){
                gradx(i, k) += NodeCo(i,j) * dphi(k, j);
            }
        }
    }
}

void GeomTriangle::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomTriangle::GetNodes(VecInt &nodes) const  {
    nodes = fNodeIndices;
}

int GeomTriangle::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomTriangle::NumNodes() {
    return nCorners;
}

GeoElementSide GeomTriangle::Neighbour(int side)  const {
    return fNeighbours[side];
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}

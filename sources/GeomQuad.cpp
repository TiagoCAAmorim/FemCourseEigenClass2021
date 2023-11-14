/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"

GeomQuad::GeomQuad() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();

    phi[0] = (1. - xi[0])*(1. - xi[1]) / 4.;
    phi[1] = (1. + xi[0])*(1. - xi[1]) / 4.;
    phi[2] = (1. + xi[0])*(1. + xi[1]) / 4.;
    phi[3] = (1. - xi[0])*(1. + xi[1]) / 4.;

    dphi(0, 0) = -(1. - xi[1]) / 4.;
    dphi(0, 1) =  (1. - xi[1]) / 4.;
    dphi(0, 2) =  (1. + xi[1]) / 4.;
    dphi(0, 3) = -(1. + xi[1]) / 4.;

    dphi(1, 0) = -(1. - xi[0]) / 4.;
    dphi(1, 1) = -(1. + xi[0]) / 4.;
    dphi(1, 2) =  (1. + xi[0]) / 4.;
    dphi(1, 3) =  (1. - xi[0]) / 4.;
}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
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

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
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

void GeomQuad::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) const {
    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    
    if (xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();
    
    phi[0] = (1. - xi[0]) / 2.;
    phi[1] = (1 + xi[0]) / 2.;
    
    dphi(0, 0) = -0.5;
    dphi(0, 1) = 0.5;
}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    
    if (xi.size() != Dimension) DebugStop(); // csi same size as dimension
    if (NodeCo.cols() != nCorners) DebugStop(); // has to supply same number of nodes as ncorners of the element
    
    // NodeCo cols: number of vertices or corners
    // NodeCo rows: number of cartesian coordinates of each corner since it could be a line living in 2d or 3d.

    int64_t nrow = NodeCo.rows();

    // Number of coordinates of the resulting x needs to be 1, 2 or 3d according to the NodeCo provided
    if (x.size() < nrow) x.resize(nrow);
    x.setZero();

    // Could have called Shape and used phi, but the mapping is quite
    // simple in this case
    for (int i = 0; i < nrow; i++) {
        x[i] = 0.5 * NodeCo(i, 0) * (1. - xi[0]) + 0.5 * NodeCo(i, 1) * (1. + xi[0]);
    }
}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    X(xi, NodeCo, x);
    
    int64_t nrow = NodeCo.rows();
    gradx.resize(nrow, 1); // nrow = number of derivatives, 1 = one point
    gradx.setZero();

    for (int i = 0; i < nrow; i++) {
        gradx(i, 0) = -0.5 * NodeCo(i, 0) + 0.5 * NodeCo(i, 1);
    }
}

void Geom1d::SetNodes(const VecInt &nodes) {
    if(nodes.rows() != 2) DebugStop();
    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt &nodes) const {
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) const {
    if(node<0 || node > 1) DebugStop();
    return fNodeIndices[node];
}

int Geom1d::NumNodes(){
    return nCorners;    
}

GeoElementSide Geom1d::Neighbour(int side) const{
    if(side <0 || side>2) DebugStop();
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    if(side < 0 || side > 2) DebugStop();
    fNeighbours[side]=neighbour;
}

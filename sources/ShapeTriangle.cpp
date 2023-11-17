//
//  ShapeTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"
#include "Shape1d.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    // Duvida: devo assumir que sempre serao 3 ou 7 lados? (lin/quad)
    DebugStop(orders.size() <= 7, "ShapeQuad::Shape: Invalid number of sides in Triangular element (>7).");
    DebugStop(orders.size() >= 3, "ShapeQuad::Shape: Invalid number of sides in Triangular element (<3).");

    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1) {
        std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    for (int i = 3; i < orders.size(); i++)
    {
        if (orders[i] > 2) {
            std::cout << "ShapeTriangle::Shape: Triangular element implemented only until order = 2. Side " << i << " has order " << orders[i] << ".\n";
            DebugStop();
        }
    }

    DebugStop(xi.size() == Dimension, "ShapeQuad::Shape: Invalid xi dimension for 2D element. Expected "+std::to_string(xi.size())+", received "+std::to_string(Dimension)+".");

    auto nf = NShapeFunctions(orders);
    phi.resize(nf);
    dphi.resize(2, nf);

    // Linear
    phi[0] =  1.-xi[0]-xi[1];
    dphi(0, 0) = -1.;
    dphi(1, 0) = -1.;

    phi[1] =  xi[0];
    dphi(0, 1) =  1.;
    dphi(1, 1) =  0.;

    phi[2] =  xi[1];
    dphi(0, 2) =  0.;
    dphi(1, 2) =  1.;

    int nfi= 3;

    // Quadratic
    if (orders.size() > 3){
        if (NShapeFunctions(3, orders[3]) > 0){
            phi[nfi] = (1.-xi[0]-xi[1]) * xi[0];
            dphi(0, nfi) = 1 - 2*xi[0] - xi[1];
            dphi(1, nfi) = -xi[0];
            nfi++;
        };
        if (NShapeFunctions(4, orders[4]) > 0){
            phi[nfi] = xi[0] * xi[1];
            dphi(0, nfi) = xi[1];
            dphi(1, nfi) = xi[0];
            nfi++;
        };
        if (NShapeFunctions(5, orders[5]) > 0){
            phi[nfi] = (1.-xi[0]-xi[1]) * xi[1];
            dphi(0, nfi) = -xi[1];
            dphi(1, nfi) = 1 - 2*xi[1] - xi[0];
            nfi++;
        };
        if (NShapeFunctions(6, orders[6]) > 0){
            phi[nfi] = (1.-xi[0]-xi[1]) * xi[0] * xi[1];
            dphi(0, nfi) = xi[1] - xi[1]*xi[1] - 2 * xi[0] * xi[1];
            dphi(1, nfi) = xi[0] - xi[0]*xi[0] - 2 * xi[0] * xi[1];
            nfi++;
        };
    }

    DebugStop(nf == nfi, "ShapeQuad::Shape: Number of shape functions does not match! Expected "+std::to_string(nf)+", got "+std::to_string(nfi)+" Check code.");
}

/// returns the number of shape functions associated with a side
int ShapeTriangle::NShapeFunctions(int side, int order){
    switch(side) {
        case 0:
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
        case 5:
            return order-1;
        case 6:
            return (order-1)*(order-1); //0;
    }

    DebugStop();
    std::cout << "ShapeTriangle::NShapeFunctions, bad parameter side " << std::endl;
    return 0;
}

/// returns the total number of shape functions
int ShapeTriangle::NShapeFunctions(VecInt &orders){

    int res=3;
    for(int in=3; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }

    return res;

}

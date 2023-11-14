//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"
#include "ShapeQuad.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){

    // Duvida: devo assumir que sempre serao 4 ou 9 lados? (lin/quad)
    DebugStop(orders.size() <= 9, "ShapeQuad::Shape: Invalid number of sides in Quad element (>9).");
    DebugStop(orders.size() >= 4, "ShapeQuad::Shape: Invalid number of sides in Quad element (<4).");

    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments (order < 0) in side" << i << ".\n";
            DebugStop();
        }
    }
        if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments (order > 1) in sides 0,1,2 or 3.\n";
        DebugStop();
    }

    for (int i = 4; i < orders.size(); i++)
    {
        if (orders[i] > 2) {
            std::cout << "ShapeQuad::Shape: Quad element implemented only until order = 2. Side " << i << " has order " << orders[i] << ".\n";
            DebugStop();
        }
    }

    DebugStop(xi.size() == Dimension, "ShapeQuad::Shape: Invalid xi dimension for 2D element. Expected "+std::to_string(xi.size())+", received "+std::to_string(Dimension)+".");

    auto nf = NShapeFunctions(orders);
    // std::cout << "ShapeQuad::NShapeFunctions, number of shape function " << nf << std::endl;
    // std::cout << "ShapeQuad::NShapeFunctions, number of elements " << orders.size() << std::endl;

    phi.resize(nf);
    dphi.resize(2, nf);

    phi[0] = (1. - xi[0])*(1. - xi[1]) / 4.;
    dphi(0, 0) = -(1. - xi[1]) / 4.;
    dphi(1, 0) = -(1. - xi[0]) / 4.;

    phi[1] = (1. + xi[0])*(1. - xi[1]) / 4.;
    dphi(0, 1) =  (1. - xi[1]) / 4.;
    dphi(1, 1) = -(1. + xi[0]) / 4.;

    phi[2] = (1. + xi[0])*(1. + xi[1]) / 4.;
    dphi(0, 2) =  (1. + xi[1]) / 4.;
    dphi(1, 2) =  (1. + xi[0]) / 4.;

    phi[3] = (1. - xi[0])*(1. + xi[1]) / 4.;
    dphi(0, 3) = -(1. + xi[1]) / 4.;
    dphi(1, 3) =  (1. - xi[0]) / 4.;

    // Duvida: Posso ter um dos lados linear e os outros quadraticos? Como fica a ordem?
    // Duvida: Como fazer com ordens superiores? Como garanto que a funcao de 3a ordem esta da direcao 'correta'?
    int nfi= 4;

    if (orders.size() > 4){
        if (NShapeFunctions(4, orders[4]) > 0){
            phi[nfi] = (1. - xi[0]*xi[0])*(1. - xi[1]) / 2.;
            dphi(0, nfi) = -xi[0] * (1. - xi[1]);
            dphi(1, nfi) = -(1. - xi[0]*xi[0]) / 2.;
            nfi++;
        };
        if (NShapeFunctions(5, orders[5]) > 0){
            phi[nfi] = (1. - xi[1]*xi[1])*(1. + xi[0]) / 2.;
            dphi(0, nfi) = (1. - xi[1]*xi[1]) / 2.;
            dphi(1, nfi) = -xi[1] * (1. + xi[0]);
            nfi++;
        };
        if (NShapeFunctions(6, orders[6]) > 0){
            phi[nfi] = (1. - xi[0]*xi[0])*(1. + xi[1]) / 2.;
            dphi(0, nfi) = -xi[0] * (1. + xi[1]);
            dphi(1, nfi) = (1. - xi[0]*xi[0]) / 2.;
            nfi++;
        };
        if (NShapeFunctions(7, orders[7]) > 0){
            phi[nfi] = (1. - xi[1]*xi[1])*(1. - xi[0]) / 2.;
            dphi(0, nfi) = -(1. - xi[1]*xi[1]) / 2.;
            dphi(1, nfi) = -xi[1] * (1. - xi[0]);
            nfi++;
        };
        if (NShapeFunctions(8, orders[8]) > 0){
            phi[nfi] = (1. - xi[0]*xi[0])*(1. - xi[1]*xi[1]);
            dphi(0, nfi) = -2*xi[0] * (1. - xi[1]*xi[1]);
            dphi(1, nfi) = -2*xi[1] * (1. - xi[0]*xi[0]);
            nfi++;
        };
    }

    DebugStop(nf == nfi, "ShapeQuad::Shape: Number of shape functions does not match! Expected "+std::to_string(nf)+", got "+std::to_string(nfi)+" Check code.");
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 7
    else if(side==8)
        return ((order-1)*(order-1));

    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();

    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){

    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }

    return res;
}

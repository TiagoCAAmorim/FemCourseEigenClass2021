

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"

using std::cout;
using std::endl;
using std::cin;

void Integrate1D();
void Integrate2DQuad();
void Integrate2DTri();

int main() {
    Integrate1D();
    Integrate2DQuad();
    Integrate2DTri();
    return 0;
}

void Integrate1D (){
    std::cout << "Teste integral 1D" << std::endl;
    auto func = [](double x){return x*x;};

    IntRule1d oned(2);
    int np = oned.NPoints();
    // oned.Print(std::cout);
    double integral = 0.;
    VecDouble co(1);
    double weight;
    for (int ip = 0; ip < np; ip++) {

        oned.Point(ip, co, weight);
        double val = func(co[0]);
        integral += val*weight;
    }
    std::cout << "  Correto = 2/3 (0.666666...),   Calculado = " << integral << std::endl;
}

void Integrate2DQuad (){
    std::cout << "Teste integral 2D: Quadrilatero" << std::endl;
    auto func = [](VecDouble xi){return xi[0]*xi[0]*xi[1]*xi[1];};

    IntRuleQuad quadrule(2);
    const int np = quadrule.NPoints();
    // quadrule.Print(std::cout);
    double integral = 0.;
    VecDouble co(2);
    double weight;
    for (int ip = 0; ip < np; ip++) {
        quadrule.Point(ip, co, weight);
        double val = func(co);
        integral += val*weight;
    }
    std::cout << "  Correto = 4/9 (0.44444...),   Calculado = " << integral << std::endl;
}

void Integrate2DTri (){
    std::cout << "Teste integral 2D: Triangulo" << std::endl;
    auto func = [](VecDouble xi){return xi[0]*xi[1]*xi[0]*xi[1];};

    IntRuleTriangle trirule(4);
    const int np = trirule.NPoints();
    // trirule.Print(std::cout);
    double integral = 0.;
    VecDouble co(2);
    double weight;
    for (int ip = 0; ip < np; ip++) {
        trirule.Point(ip, co, weight);
        double val = func(co);
        // printf("  %d. f(%.4f,%.4f) = %.5f, w = %.4f\n",ip, co[0], co[1], val, weight);
        integral += val*weight;
    }
    std::cout << "  Correto = 1/180 (0.0055555...),   Calculado = " << integral << std::endl;
}
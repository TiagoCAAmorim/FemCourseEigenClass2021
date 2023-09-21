

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

int main() {
    Integrate1D();
    return 0;
}

void Integrate1D (){

    // test an integration rule
    // lambda expression
    auto func = [](double x){return x*x;};
    
    IntRule1d oned(2);
    int np = oned.NPoints();
    double integral = 0.;
    VecDouble co(1);
    double weight;
    for (int ip = 0; ip < np; ip++) {
        
        oned.Point(ip, co, weight);
        double val = func(co[0]);
        integral += val*weight;
    }
    std::cout << "espera se 2/3 obtem se " << integral << std::endl;    
}

#include <iostream>
#include <math.h>
#include "VTKGeoMesh.h"
#include "Geom1d.h"

using namespace std;


void ShapeFunc1D();

int main() {
    ShapeFunc1D();
    return 0;
}

void ShapeFunc1D(){
    
    Geom1d geom1d;
    
    // Number of nodes of the element and dimension where the element exists
    const int64_t nnodes = 2;
    const int dim = 1;
            
    // Element nodes
    MatrixDouble NodeCo(dim,nnodes);
    NodeCo(0,0) = 0.;
    NodeCo(0,1) = 4.;
    
    // Testing for several points
    VecDouble xi(dim), x(dim);
    MatrixDouble gradx(dim,1);
    VecDouble testpts(5);
    testpts << -1,-0.5,0.,0.5,1.;
    for(int i = 0 ; i < testpts.size() ; i++) {
        xi[0] = testpts[i];
        geom1d.GradX(xi, NodeCo, x, gradx);
        cout << "------ csi = " << xi[0] << " ------" << endl;
        cout << "x = " << x << "\t|\tgradx = " << gradx << endl;
    }    
}

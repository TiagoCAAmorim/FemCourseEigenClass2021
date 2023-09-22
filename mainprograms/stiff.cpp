#include <iostream>
#include <math.h>
#include <string>
#include "VTKGeoMesh.h"
#include "Geom1d.h"
#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "CompMesh.h"
#include "GeoElementTemplate.h"
#include "GeoMesh.h"
#include "Poisson.h"

using namespace std;


void CallStiffFunc();

int main() {
    CallStiffFunc();
    return 0;
}

void CallStiffFunc(){

    std::cout << "'Horizontal beam'\n";
    //Nodes in the geometry mesh
    const int dim = 2;
    const int numnodes = 2;
    GeoMesh gmesh;
    gmesh.SetDimension(dim);
    gmesh.SetNumNodes(numnodes);
    GeoNode gnod0, gnod1;
    VecDouble co0(dim), co1(dim);
    co0 << 0., 0.;
    co1 << 4., 0.;
    gnod0.SetCo(co0);
    gnod1.SetCo(co1);
    gmesh.Node(0) = gnod0;
    gmesh.Node(1) = gnod1;
    
    //New 1D element between nodes 0 and 1
    int materialid = 0;
    VecInt nodeindices(2);
    nodeindices << 0,1;
    int index = 0;
    GeoElementTemplate<Geom1d> geo(nodeindices,materialid,&gmesh,index);

    //New computational mesh based on geometry mesh, with linear functions (?)
    int order = 1;
    CompMesh cmesh(&gmesh);
    cmesh.SetDefaultOrder(order);
    
    // Criando material
    MatrixDouble perm(3,3);
    perm.setIdentity();
    perm(2,2) = 2.;
    Poisson poi(materialid, perm);
    cmesh.SetMathStatement(materialid, &poi);
    
    CompElementTemplate<Shape1d> cel(index,&cmesh,&geo);
    MatrixDouble ek(order+1,order+1),ef(order+1,1);
    cel.CalcStiff(ek, ef);
    cout << ek << endl;
}
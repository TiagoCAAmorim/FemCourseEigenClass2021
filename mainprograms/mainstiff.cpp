#include <iostream>
#include <math.h>
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
    
    
    const int dim = 1;
    const int numnodes = 2;
    GeoMesh gmesh;
    gmesh.SetDimension(dim);
    gmesh.SetNumNodes(numnodes);
    GeoNode gnod0, gnod1;
    VecDouble co0(1), co1(1);
    co0 << 0.;
    co1 << 4.;
    gnod0.SetCo(co0);
    gnod1.SetCo(co1);
    
    gmesh.Node(0) = gnod0;
    gmesh.Node(1) = gnod1;
    
    int materialid = 0;
    VecInt nodeindices(2);
    nodeindices << 0,1;
    int index = 0;
//    const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index
    GeoElementTemplate<Geom1d> geo(nodeindices,materialid,&gmesh,index);
        
    CompMesh cmesh(&gmesh);
    cmesh.SetDefaultOrder(1);
    
    // Criando material
    MatrixDouble perm(3,3);
    perm.setIdentity();
//    perm(1,1) = 2.;
    Poisson poi(materialid, perm);
    cmesh.SetMathStatement(materialid, &poi);
    
    CompElementTemplate<Shape1d> cel(index,&cmesh,&geo);
    MatrixDouble ek(2,2),ef(2,1);
    cel.CalcStiff(ek, ef);
    cout << ek << endl;
}
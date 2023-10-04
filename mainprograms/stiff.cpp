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

    std::cout << "Setting up geometric mesh\n";
    const int dim = 2;
    const int numnodes = 3;
    GeoMesh gmesh;
    gmesh.SetDimension(dim);
    gmesh.SetNumNodes(numnodes);
    GeoNode gnod0, gnod1, gnod2;
    VecDouble co0(dim), co1(dim), co2(dim);
    co0 << 0., 0.;
    co1 << 4., 0.;
    co2 << 0., 4.;
    gnod0.SetCo(co0);
    gnod1.SetCo(co1);
    gnod2.SetCo(co2);
    gmesh.Node(0) = gnod0;
    gmesh.Node(1) = gnod1;
    gmesh.Node(2) = gnod2;
    
    std::cout << "Setting up material\n";
    int materialid = 0;
    MatrixDouble perm(3,3);
    perm.setIdentity();
    perm(1,1) = 2.;
    Poisson poi(materialid, perm);

    std::cout << "Setting up 'Horizontal beam': #0\n";
    VecInt nodeindices_h(2);
    nodeindices_h << 0,1;
    int index_h = 0;
    // Builds element and adds to geometric mesh
    GeoElementTemplate<Geom1d> geo_h(nodeindices_h,materialid,&gmesh,index_h);

    std::cout << "Setting up 'Vertical beam': #1\n";
    VecInt nodeindices_v(2);
    nodeindices_v << 0,2;
    int index_v = 1;
    GeoElementTemplate<Geom1d> geo_v(nodeindices_v,materialid,&gmesh,index_v);

    std::cout << "Building computational mesh: linear\n";
    int order = 1;
    CompMesh cmesh(&gmesh);
    cmesh.SetDefaultOrder(order);
    cmesh.SetMathStatement(materialid, &poi); //Adds material to computational mesh
    
    MatrixDouble ek(order+1,order+1),ef(order+1,1);

    std::cout << "Stiffness matrix for 'Horizontal beam'\n";
    CompElementTemplate<Shape1d> cel_h(index_h,&cmesh,&geo_h); // Creates computational element and adds to Computational mesh
    cel_h.CalcStiff(ek, ef);
    cout << ek << endl;

    std::cout << "Stiffness matrix for 'Vertical beam'\n";
    CompElementTemplate<Shape1d> cel_v(index_v,&cmesh,&geo_v);
    cel_v.CalcStiff(ek, ef);
    cout << ek << endl;


    std::cout << "Building computational mesh: quadratic\n";
    order = 2;
    CompMesh cmesh_quad(&gmesh);
    cmesh_quad.SetDefaultOrder(order);
    cmesh_quad.SetMathStatement(materialid, &poi);
    
    MatrixDouble ek_quad(order+1,order+1),ef_quad(order+1,1);

    std::cout << "Stiffness matrix for 'Horizontal beam'\n";
    CompElementTemplate<Shape1d> cel_h_quad(index_h,&cmesh_quad,&geo_h);
    cel_h_quad.CalcStiff(ek_quad, ef_quad);
    cout << ek_quad << endl;

    std::cout << "Stiffness matrix for 'Vertical beam'\n";
    CompElementTemplate<Shape1d> cel_v_quad(index_v,&cmesh_quad,&geo_v);
    cel_v_quad.CalcStiff(ek_quad, ef_quad);
    cout << ek_quad << endl;
}
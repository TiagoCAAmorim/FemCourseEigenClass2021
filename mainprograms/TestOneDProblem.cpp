//
//  TestOneDProblem.cpp
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
        Os testes foram preparados com um proposito educacional,
        recomenda-se que o aluno entenda a funcionalidade de cada
        teste e posteriormente use com seu código caso a caso
*/
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
#include <iostream>
#include <math.h>
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "MathStatement.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "IntRule.h"
#include "PostProcessTemplate.h"
#include "Poisson.h"

using std::cout;
using std::endl;
using std::cin;

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);

int main ()
{
    // Geometric mesh created, with data read from a Gmsh file format: nodes, elements, material id
    GeoMesh gmesh;
    ReadGmsh read;
    std::string filename("D:/FemCourseEigenClass2021/mainprograms/oneD.msh");
#ifdef MACOSX
    filename = "../"+filename;
#endif
    read.Read(gmesh,filename);

    // Created comp mesh and defined mathstatement 
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    mat1->SetDimension(1); // This mathstatement will be applied to 1D elements
    
    // Lambda function to define right-side of differential equation
    auto force = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 1.; // Poisson: -u'' = 1; 'Standard problem': res[0] = x[0]
    };
    mat1->SetForceFunction(force);

    // Boundary conditions
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    // L2 projection solves u=f => find u that approximates F in the function space
    // 0 => Dirichilet 
    L2Projection *bc_linha = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,3,proj,val1,val2);

    // Defining mathstatements as a vector and passing to computational mesh
    std::vector<MathStatement *> mathvec = {nullptr,mat1,bc_linha,bc_point};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(1);

    // Creates the computational elements
    cmesh.AutoBuild();
    cmesh.Resequence();
   
    // Creates matrices by element, global matrices and solves problem
    Analysis AnalysisLoc(&cmesh);
    AnalysisLoc.RunSimulation();
    
    // Post-porcessing: comparison to exact result
    PostProcessTemplate<Poisson> postprocess;
    postprocess.SetExact(exact);
    
    VecDouble errvec;
    errvec = AnalysisLoc.PostProcessError(std::cout, postprocess);
        
    return 0;
}
void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv){

    deriv(0,0) = 4-point[0];
    val[0]=point[0]*(8.-point[0])/2.;
    return;
}



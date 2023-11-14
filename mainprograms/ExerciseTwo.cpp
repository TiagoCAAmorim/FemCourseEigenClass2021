//  ExerciseTwo.cpp
//      Based on TestTwoDProblem.cpp
//      Created by Tiago Amorim, 13/11/2023.
/*
    Exercise 2 in FEM course
    Know responses: 1, x, y, xy, x2, y2, x2y2
    Test each with four meshes:
        1 quadrilateral
        4 quadrilaterals
        2 triangles
        8 triangles
    Each test with linear and quadratic elements.
    7 x 4 x 2 = 56 tests!!
    Domain: [{0,0},{1,1}]
 */
//
#include <iostream>
#include <math.h>
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "PostProcessTemplate.h"

void f1_1quad_lin(){
    GeoMesh gmesh;
    ReadGmsh read;
    std::string filename("D:/FemCourseEigenClass2021/ex02/Mesh1Quad.msh");
#ifdef MACOSX
    filename = "../"+filename;
#endif
    read.Read(gmesh,filename);

    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    mat1->SetDimension(1);

    auto force = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 0;
    };
    mat1->SetForceFunction(force);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,3,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_point,bc_linha};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(1);
    cmesh.AutoBuild();
    cmesh.Resequence();

    Analysis locAnalysis(&cmesh);
    locAnalysis.RunSimulation();
    PostProcessTemplate<Poisson> postprocess;
    auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = (1.-x[0])*x[0]*(1-x[1])*x[1];
        deriv(0,0) = (1.-2.*x[0])*(1-x[1])*x[1];
        deriv(1,0) = (1-2.*x[1])*(1-x[0])*x[0];
    };

//    if (!strcmp("Sol", name.c_str())) return ESol;
//    if (!strcmp("DSol", name.c_str())) return EDSol;
//    if (!strcmp("Flux", name.c_str())) return EFlux;
//    if (!strcmp("Force", name.c_str())) return EForce;
//    if (!strcmp("SolExact", name.c_str())) return ESolExact;
//    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;
    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");
    postprocess.AppendVariable("Flux");
    postprocess.AppendVariable("Force");
    postprocess.AppendVariable("SolExact");
    postprocess.AppendVariable("DSolExact");
    postprocess.SetExact(exact);
    mat1->SetExactSolution(exact);
    locAnalysis.PostProcessSolution("D:/FemCourseEigenClass2021/ex02/f1_1quad_lin.vtk", postprocess);
    // locAnalysis.PostProcessSolution("quads.vtk", postprocess);

    VecDouble errvec;
    errvec = locAnalysis.PostProcessError(std::cout, postprocess);
}


int main(){
    f1_1quad_lin();
    return 0;
}
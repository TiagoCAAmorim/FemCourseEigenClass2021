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

void test2D(std::string filename, std::string functionname, const std::function<void(const VecDouble &co, VecDouble &result)> &force, const std::function<void(const VecDouble &loc, VecDouble &result, MatrixDouble &deriv)> &exact, int order){
    GeoMesh gmesh;
    ReadGmsh read;
    read.Read(gmesh,filename+".msh");

    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    mat1->SetDimension(2);
    mat1->SetForceFunction(force);

    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_left   = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_right  = new L2Projection(0,3,proj,val1,val2);
    L2Projection *bc_bottom = new L2Projection(0,4,proj,val1,val2);
    L2Projection *bc_top    = new L2Projection(0,5,proj,val1,val2);
    bc_left->SetExactSolution(exact);
    bc_right->SetExactSolution(exact);
    bc_bottom->SetExactSolution(exact);
    bc_top->SetExactSolution(exact);

    std::vector<MathStatement *> mathvec = {0,mat1,bc_left, bc_right, bc_bottom, bc_top};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(order);
    cmesh.AutoBuild();  // Monta os DOFs e as equacoes associadas.
    cmesh.Resequence();  //Duvida: ta' fazendo resequence duas vezes. Ja' tem dentro de AutoBuild.

    Analysis locAnalysis(&cmesh);
    locAnalysis.RunSimulation();
    PostProcessTemplate<Poisson> postprocess;

    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");
    postprocess.AppendVariable("Flux");
    postprocess.AppendVariable("Force");
    postprocess.AppendVariable("SolExact");
    postprocess.AppendVariable("DSolExact");
    postprocess.SetExact(exact);
    mat1->SetExactSolution(exact);
    locAnalysis.PostProcessSolution(filename+"_"+functionname+"_"+std::to_string(order)+".vtk", postprocess);

    VecDouble errvec;
    errvec = locAnalysis.PostProcessError(std::cout, postprocess);
}

void mutipletests(std::string functionname, const std::function<void(const VecDouble &co, VecDouble &result)> &force, const std::function<void(const VecDouble &loc, VecDouble &result, MatrixDouble &deriv)> &exact){
    std::string foldername = "D:/FemCourseEigenClass2021/ex02/";
    std::vector<std::string> filenames = {"Mesh1Quad", "Mesh4Quad", "Mesh2Tri", "Mesh8Tri"};
    int max_order = 2; // (1 or 2)

    for (int i=0; i<filenames.size(); i++){
        for (int j=0; j<max_order; j++){
            std::cout << "" << std::endl;
            std::cout << "##### " << filenames[i] << " #####" << std::endl;
            std::cout << "#####   exact function: " << functionname << " #####" << std::endl;
            std::cout << "#####   order: " << (j+1) << " #####" << std::endl;
            test2D(foldername+filenames[i], functionname, force, exact, j+1);
        }
    }
}

void mutipletests2(std::string functionname, const std::function<void(const VecDouble &co, VecDouble &result)> &force, const std::function<void(const VecDouble &loc, VecDouble &result, MatrixDouble &deriv)> &exact){
    std::string foldername = "D:/FemCourseEigenClass2021/ex02/";
    std::vector<std::string> filenames = {"Mesh16Quad", "Mesh64Quad", "Mesh256Quad", "Mesh32Tri", "Mesh128Tri", "Mesh512Tri"};
    int max_order = 1; // (1 or 2)

    for (int i=0; i<filenames.size(); i++){
        for (int j=0; j<max_order; j++){
            std::cout << "" << std::endl;
            std::cout << "##### " << filenames[i] << " #####" << std::endl;
            std::cout << "#####   exact function: " << functionname << " #####" << std::endl;
            std::cout << "#####   order: " << (j+1) << " #####" << std::endl;
            test2D(foldername+filenames[i], functionname, force, exact, j+1);
        }
    }
}


int main(){
    SetDebug(false);

    auto exact1 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = 1;
        deriv(0,0) = 0;
        deriv(1,0) = 0;
    };
    auto force1 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 0;
    };
    mutipletests("unity", force1, exact1);

    auto exact2 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0];
        deriv(0,0) = 1;
        deriv(1,0) = 0;
    };
    auto force2 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 0;
    };
    mutipletests("x", force2, exact2);

    auto exact3 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[1];
        deriv(0,0) = 0;
        deriv(1,0) = 1;
    };
    auto force3 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 0;
    };
    mutipletests("y", force3, exact3);

    auto exact4 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0]*x[1];
        deriv(0,0) = x[1];
        deriv(1,0) = x[0];
    };
    auto force4 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 0;
    };
    mutipletests("xy", force4, exact4);

    auto exact5 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0]*x[0];
        deriv(0,0) = 2*x[0];
        deriv(1,0) = 0.;
    };
    auto force5 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = -2;
    };
    mutipletests("x2", force5, exact5);

    auto exact6 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[1]*x[1];
        deriv(0,0) = 0.;
        deriv(1,0) = 2*x[1];
    };
    auto force6 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = -2;
    };
    mutipletests("y2", force6, exact6);

    auto exact7 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0]*x[0]*x[1]*x[1];
        deriv(0,0) = 2*x[0]*x[1]*x[1];
        deriv(1,0) = 2*x[0]*x[0]*x[1];
    };
    auto force7 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = -2*x[0]*x[0]-2*x[1]*x[1];
    };
    mutipletests("x2y2", force7, exact7);
    // mutipletests2("x2y2", force7, exact7);

    return 0;
}
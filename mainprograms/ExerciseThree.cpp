//  ExerciseTwo.cpp
//      Based on TestTwoDProblem.cpp
//      Created by Tiago Amorim, 13/11/2023.
/*
    Exercise 3 in FEM course
    Known response:
    Test each with 10 meshes:
        quadrilaterals x 5 refinements
        triangles x 5 refinements
    Each test with linear and quadratic elements.
    10 x 2 = 20 tests
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
    L2Projection *bc   = new L2Projection(0,2,proj,val1,val2);
    bc->SetExactSolution(exact);

    std::vector<MathStatement *> mathvec = {0,mat1,bc};
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
    std::string foldername = "D:/FemCourseEigenClass2021/ex03/";
    std::vector<std::string> filenames = {"MalhaAnelQuad_0", "MalhaAnelTria_0"};
    int max_refine = 5; // (1 to 5)
    int max_order = 2; // (1 or 2)

    for (int i=0; i<filenames.size(); i++){
        for (int j=0; j<max_order; j++){
            for (int k=0; k<max_refine; k++){
                std::cout << "" << std::endl;
                std::cout << "##### " << filenames[i] << (k+1) << " #####" << std::endl;
                std::cout << "#####   exact function: " << functionname << " #####" << std::endl;
                std::cout << "#####   order: " << (j+1) << " #####" << std::endl;
                test2D(foldername+filenames[i]+std::to_string(k+1), functionname, force, exact, j+1);
            }
        }
    }
}

int main(){
    SetDebug(false);

    auto exactX = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0];
        deriv(0,0) = 1.;
        deriv(1,0) = 0.;
    };
    auto forceX = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 0;
    };
    // mutipletests("x", forceX, exactX);

    auto exactXY = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0]*x[1];
        deriv(0,0) = x[1];
        deriv(1,0) = x[0];
    };
    auto forceXY = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 0;
    };
    // mutipletests("xy", forceXY, exactXY);

    auto exactX2 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0]*x[0];
        deriv(0,0) = 2*x[0];
        deriv(1,0) = 0.;
    };
    auto forceX2 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = -2;
    };
    // mutipletests("x2", forceX2, exactX2);

    auto exactX2Y2 = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0]*x[0]*x[1]*x[1];
        deriv(0,0) = 2*x[0]*x[1]*x[1];
        deriv(1,0) = 2*x[0]*x[0]*x[1];
    };
    auto forceX2Y2 = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = -2*x[0]*x[0]-2*x[1]*x[1];
    };
    // mutipletests("x2y2", forceX2Y2, exactX2Y2);

    auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        double pi = 3.1415926535897932384626433832795028842;
        val[0] = sin(3. * pi * x[0]) * sin(pi * x[1]);
        deriv(0,0) = 3. * pi * cos(3. * pi * x[0]) * sin(pi * x[1]);
        deriv(1,0) = pi * sin(3. * pi * x[0]) * cos(pi * x[1]);
    };
    auto force = [](const VecDouble &x, VecDouble &res)
    {
        double pi = 3.1415926535897932384626433832795028842;
        res[0] = 10. * pi * pi * sin(3. * pi * x[0]) * sin(pi * x[1]);
    };
    mutipletests("Ex03", force, exact);

    return 0;
}
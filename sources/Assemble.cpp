/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Assemble.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "MathStatement.h"
#include "CompElement.h"



Assemble::Assemble() {
}

Assemble::Assemble(CompMesh *mesh) {
    cmesh = mesh;
}

Assemble::Assemble(const Assemble &copy) {
    cmesh = copy.cmesh;
}

Assemble &Assemble::operator=(const Assemble &copy) {
    cmesh = copy.cmesh;
    return *this;
}

void Assemble::SetMesh(CompMesh *mesh) {
    cmesh = mesh;
}

int64_t Assemble::NEquations() {
    int64_t neq = 0;
    int64_t i, ncon = cmesh->GetDOFVec().size();
    for (i = 0; i < ncon; i++) {
        DOF dof = cmesh->GetDOF(i);
        int64_t dofsize = dof.GetNShape() * dof.GetNState();
        neq += dofsize;
    }
    return neq;
}

void Assemble::OptimizeBandwidth() {
}

void Assemble::Compute(SparseMat &globmat, MatrixDouble &rhs) {

    auto neq = NEquations();

    globmat.resize(neq, neq); // Duvida: nao
    globmat.setZero();
    rhs.resize(neq, 1);
    rhs.setZero();

    int64_t nelem = cmesh->GetGeoMesh()->NumElements();
    for (int el = 0; el < nelem; el++) {
        CompElement *cel = cmesh->GetElement(el);

        int nshape = cel->NShapeFunctions();
        int nstate = cel->GetStatement()->NState(); // Number of state variables (right-side function dimension)
        int neqs = nstate * nshape;
        MatrixDouble ek(neqs, neqs);
        MatrixDouble ef(neqs, 1);
        ek.setZero();
        ef.setZero();

        std::cout << "    Element #" << (el+1) << ": " << nshape << " shape function(s)." << std::endl;
        cel->CalcStiff(ek, ef);

        VecInt destIndex(neqs);
        int count = 0;
        destIndex.setZero();
        const int64_t ndof = cel->NDOF();
        for (int i=0; i<ndof; i++){
            DOF& dofi = cel->GetDOF(i);
            const int64_t firstEq = dofi.GetFirstEquation();
            const int64_t neqDofi = dofi.GetNShape() * dofi.GetNState();
            for (int idofi=0; idofi<neqDofi; idofi++){
                destIndex[count++] = firstEq + idofi;
            }
        }

        for (int i=0; i<neqs; i++){
            int64_t dest_i = destIndex[i];
            rhs(dest_i,0) += ef(i,0);
            for (int j=0; j<neqs; j++){
                int64_t dest_j = destIndex[j];
                globmat.coeffRef(dest_i,dest_j) += ek(i,j);
            }
        }


    }
}

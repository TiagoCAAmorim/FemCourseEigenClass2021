/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream>
///\endcond
#include "IntRule1d.h"
#include "IntRuleQuad.h"

IntRuleQuad::IntRuleQuad(){
    SetOrder(0);
}

IntRuleQuad::IntRuleQuad(int order) {
    SetOrder(order);
}

void IntRuleQuad::SetOrder(int order) {
    fOrder = order;
    const int np1d = order/2 + 1;
    gaulegQuad(-1,1,fPoints, fWeights, np1d);
}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, MatrixDouble &co, VecDouble &w, const int np1d) {
    IntRule1d x;
    IntRule1d y;

    int n = np1d;

    VecDouble cox(n);
    VecDouble coy(n);
    VecDouble wx(n);
    VecDouble wy(n);

    x.gauleg(x1, x2, cox, wx);
    y.gauleg(x1, x2, coy, wy);

    co.resize(n*n, 2);
    w.resize(n * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            co(j+i*n, 0) = cox[i];
            co(j+i*n, 1) = coy[j];
            w[j+i*n] = wx[i] * wy[j];
        }
    }
}

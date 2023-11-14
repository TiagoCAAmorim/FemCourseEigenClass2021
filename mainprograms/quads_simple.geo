//+
SetFactory("OpenCASCADE");

elements = 1; // numero de elementos (alterar para fazer o projeto)
nodes = elements + 1; // numero de nós (gerando malhas com mesmo numero de nos na base e na altura)
recombine = 1; // fator para decidir se a malha é triangular (0) ou quadricular (1)

Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Point("fix", 3) = {1};
//+
Physical Curve("contorno", 2) = {4, 1, 2, 3};
//+
Physical Surface("plano", 1) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 2, 3, 4} = nodes Using Progression 1;
//Transfinite Curve {4, 2} = 2 Using Progression 1;
//+
//Transfinite Curve {1, 3} = 2 Using Progression 1;

If (recombine)
	Recombine Surface {1};
EndIf
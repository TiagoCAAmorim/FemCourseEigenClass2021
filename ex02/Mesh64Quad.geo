// ++++++++++++++++++++++++++++++ //
// 	DADOS DE ENTARDA DO CODIGO //
// ++++++++++++++++++++++++++++++ //

// Kernel do GMSH
SetFactory("OpenCASCADE");

// Variaveis de entrada do codigo geo
elements = 8; // numero de elementos (alterar para fazer o projeto)
nodes = elements + 1; // numero de nós (gerando malhas com mesmo numero de nos na base e na altura)
recombine = 1; // fator para decidir se a malha é triangular (0) ou quadricular (1)
lc = 5e-1;

// ++++++++++++++++++++++++++++++ //
// 	DEFININDO A GEOMETRIA DO DOMINIO //
// ++++++++++++++++++++++++++++++ //
// 			NAO ALTERAR

// Pontos
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// Linhas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Curvas e Superficies
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

// ++++++++++++++++++++++++++++++ //
// 	DEFININDO AS CONDICOES DE CONTORNO //
// ++++++++++++++++++++++++++++++ //

// Organizadas em "nome", "materialID"
Physical Surface("Domain", 1) = {1};

Physical Curve("BCleft", 2) = {4};
Physical Curve("BCright", 3) = {2};
Physical Curve("BCbottom", 4) = {1};
Physical Curve("BCtop", 5) = {3};

// ++++++++++++++++++++++++++++++ //
// 	DEFININDO OS PARAMETROS DA MALHA //
// ++++++++++++++++++++++++++++++ //
// 			NAO ALTERAR

// Gerando malha estruturada
Transfinite Curve {1, 2, 3, 4} = nodes Using Progression 1;
Transfinite Surface {1};

// Definindo tipo de elementos
// se recombine, então elementos quadriculares
// caso contrario, elementos triangulares (default do gmsh)
If (recombine)
	Recombine Surface {1};
EndIf
//

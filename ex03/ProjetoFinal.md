# Projeto Final

Nesta tarefa, deverão utilizar o código C++ que desenvolveram para analisar problemas bidimensionais e tirar taxas de convergência no caso de uma solução manufaturada fornecida e um domínio não trivial.

### Definição do problema

* Domínio de um quarto de anel definido com arquivo .geo fornecido.
* Equação de Poisson.
* Solução manufaturada fornecida.
* Condição de contorno de Dirichlet em todo o contorno.

### O que é pedido

* 4 conjuntos de simulações em que cada conjunto usa os seguintes elementos: quadriláteros lineares, quadriláteros quadráticos, triângulos lineares e triângulos quadráticos.
* Geração de 5 malhas com refinamentos uniformes subsequentes usando primeiro quadriláteros, e depois triângulos (total de 10 malhas). O software Gmsh deve fazer isso automaticamente utilizando o commando “Refine by splitting”.
* Como simulações devem ser rodadas utilizando elementos de ordem linear (p=1) e quadrática (p=2), apesar de gerarem 10 malhas, haverão um total de 20 simulações.
* Implementação de termo fonte e adicionar o mesmo no _MathStatement_ de Poisson.
* Identificar no código onde suprir a solução exata para cálculo do erro e chamar a função que calcula os erros. Notem que a função que calcula o erro deve calcular os mesmos nas normas L2 de $u$, L2 do gradiente de $u$, e norma H1. Portanto cada simulação deve soltar 3 números que são esses erros nas respectivas normas.
* Solução exata a ser aproximada: $u(x,y) = \sin(3 \pi x) \times \sin(\pi y)$.
* Chamar função que gera arquivo de saída com solução em formato VTK e visualização da solução aproximada no Paraview utilizando o filtro _warp by scalar_.
* Plotagem das curvas de convergência para cada conjunto de simulações. Serão um total de 4 curvas para cada norma de erro. Pede-se que separem em três gráficos, onde cada gráfico é relativo somente a uma norma de erro.

### No relatório deve constar

* Definição matemática adequada do problema analisado.
* Figuras das malhas utilizadas.
* Figuras da solução aproximada utilizando o filtro _warp by scalar_ do Paraview. No filtro _warp by scalar_, use sempre o mesmo fator de escala e indique no relatório qual valor foi usado.
* Gráficos com as curvas de convergência (erro x h), identificando quais as taxas obtidas entre as duas últimas malhas de cada uma das curvas. Note que serão 3 gráficos (relativos às normas), com 4 curvas cada um (relativas aos tipos de elemento utilizados).
* Conclusões se as taxas obtidas são as esperadas pela teoria e o porquê.
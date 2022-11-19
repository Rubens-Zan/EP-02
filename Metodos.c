#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int encontraMax(SistLinear_t *SL, int i)
{
  real_t max = SL->A[i][i]; // guarda o maior valor como o primeiro elemento da linha
  int max_index = i;        // guarda o índice do maior

  // loop para percorrer as linhas da matriz
  for (int l = i; l < SL->n; l++)
  {
    // Verifica se o elemento encontrado é maior que o maior até então
    if (ABS(SL->A[l][i]) > max)
    {
      max = SL->A[l][i];
      max_index = l;
    }
  }
  return max_index;
}

void trocaLinha(SistLinear_t *SL, int i, int max_index)

{
  // Percorre todos os elementos da linha, mudando as colunas 'k'
  for (int k = 0; k <= SL->n; k++)
  {
    real_t temp = SL->A[i][k];
    SL->A[i][k] = SL->A[max_index][k];
    SL->A[max_index][k] = temp;
  }

  // Troca o termo independente do vetor de termos independentes
  real_t temp = SL->b[max_index];
  SL->b[max_index] = SL->b[i];
  SL->b[i] = temp;
}

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal)
{
  /* para cada linha a partir da primeira */
  for (int i = 0; i < SL->n; ++i)
  {
    unsigned int max_index = encontraMax(SL, i);
    if (i != max_index)
      trocaLinha(SL, i, max_index);

    for (int k = i + 1; k < SL->n; ++k)
    {
      real_t m = SL->A[k][i] / SL->A[i][i]; // pivô

      // zera o elemento para diminuir o número de operações
      SL->A[k][i] = 0.0;

      for (int j = i + 1; j < SL->n; ++j)
        SL->A[k][j] -= SL->A[i][j] * m; //
      x[k] -= x[i] * m;
    }

    // retrosubstituição
    for (int i = SL->n - 1; i >= 0; i--)
    {
      real_t sum = 0;
      for (int j = i; j < SL->n; j++)
      {
        sum = sum + SL->A[i][j] * SL->b[j];
      }
      x[i] = (SL->A[i][SL->n + 1] - sum) / SL->A[i][i];
    }
  }
}

// VERIFICAR!!
real_t calcularNormaMaxErroAbsoluto(real_t *aux, real_t *x, int tam)
{
  real_t max = ABS(x[0] - aux[0]);

  // Percorre todos os termos
  for (int i = 1; i < tam; ++i)
  {
    // Se o módulo da subtração é maior que max
    if (ABS(x[i] - aux[i]) > max)
    {
      max = ABS(x[i] - aux[i]);
    }
  }
  return max;
}
// VERIFICAR!!


/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
*/
int gaussSeidel(SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  real_t x_calculated, error_calculated;
  int sum, num_iterations;

  // For que preenche o vetor solução com o chute inicial, "0.00"
  for (int i = 0; i < SL->n; i++)
  {
    x[i] = 0.00;
  }

  // gauss-Seidel method
  // for (num_iterations = 1; num_iterations <= MAXIT; num_iterations++)
  // {
    for (int i = 0; i < SL->n - 1 && num_iterations <= MAXIT; i++)
    {
      sum = 0;
      for (int j = 0; j < SL->n - 1; j++)
      {
        if (j != i)
        {
          // linha que calcula a soma dos coeficientes multiplicado pelo chute inicial
          sum += SL->A[i][j] * x[j];
        }
      }

      // linha que calcula o valor de c(n) - sum / a(n,n)
      x_calculated = (SL->b[i] - sum) / SL->A[i][i];

      // calcula erro em relação ao vetor de solução
      error_calculated = ABS(x[i] - x_calculated);
   printf ("erro calco : %10g    %d\n", error_calculated,num_iterations);
      if (error_calculated > erro)
      {
        x[i] = x_calculated;
      }
      else
      {
        printf("Number of iterations: %d", num_iterations);
        return 0;
      }

      num_iterations++;
    }
  // }
}

/*!
  \brief Essa função calcula a multiplicação de A*X da matriz

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param solution Solução

*/
void multiMatrix(SistLinear_t *SL, real_t *x, real_t *solution)
{
  // solution = [A]*[X]
  for (int i = 0; i < SL->n; i++)
    solution[i] = 0;

  for (int i = 0; i < SL->n; ++i)
  {
    for (int j = 0; j < SL->n; ++j)
      solution[i] += SL->A[i][j] * x[j];
  }
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param residue Valor do resíduo

  \return Norma L2 do resíduo

*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *residue)
{
  real_t norma = 0;

  // A*X
  real_t *matrix_solution;
  // alocação de memória para variável Auxiliar que guarda o valor da multiplicação da matriz
  matrix_solution = (real_t *)calloc(SL->n, sizeof(real_t));

  // função que calcula a Multiplicação A*X da matriz
  multiMatrix(SL, x, matrix_solution);

  // r0 = b0 - x0, r1 = b1 - x1, r2 = b2 - x2....
  for (int i = 0; i < SL->n; i++)
    residue[i] = SL->b[i] - matrix_solution[i];

  // Calcula da norma -> norma = sqrt(r1²+r2²+r3²+...)
  for (int i = 0; i < SL->n; i++)
    norma += residue[i] * residue[i];
  norma = sqrt(norma);
  return norma;
}

/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int refinamento(SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  // inicializa tempo da função
  double time = timestamp();

  // inicializa contador de iterações
  int iteracao = 0;

  // inicializa variáveis necessárias para que seja realizado refinamento
  real_t *aux, *w, *r, norma, *residue;

  // alocação de memória para variável auxiliar "aux"
  aux = (real_t *)calloc(SL->n, sizeof(real_t));

  // alocação de memória para variável que guarda valor do resíduo "residue"
  residue = (real_t *)calloc(SL->n, sizeof(real_t));

  // alocação de memória para variável "w" erro
  w = (real_t *)calloc(SL->n, sizeof(real_t));

  // alocação de memória para variável "resíduo" atual
  r = (real_t *)calloc(SL->n, sizeof(real_t));

  // guarda valor da Norma calculado pela função
  norma = normaL2Residuo(SL, x, residue);

  while (norma > 5.0 && iteracao < MAXIT)
  {
    iteracao++;

    multiMatrix(SL, x, aux);

    for (int i = 0; i < SL->n; i++)
      r[i] = SL->b[i] - aux[i];

    // copia sistema linear original para nova iteração do refinamento
    SistLinear_t *newSL = alocaSisLin(SL->n);

    // guarda SL antigo em nova variável para atualização do resíduo
    for (int i = 0; i < SL->n; i++)
      for (int j = 0; j < SL->n; j++)
        newSL->A[i][j] = SL->A[i][j];

    // vetor de termos indepentes que recebe valor do resíduo[i]
    for (int i = 0; i < SL->n; i++)
      newSL->b[i] = r[i];

    // ajusta variável tamanho
    newSL->n = SL->n;

    /*
    A PENSAR... newSL->erro = SL->erro;
    */

    // Chamada da eliminação de Gauss que resolve-> A*w = r
    eliminacaoGauss(newSL, w, tTotal);

    for (int i = 0; i < SL->n; i++)
    {
      x[i] += w[i];
      w[i] = 0;
    }

    norma = normaL2Residuo(SL, x, residue);
  }

  // Coloca tempo decorrido para execução da função no vetor tTotal do programa
  time = timestamp() - time;
  tTotal[0] = time;

  return iteracao;
}

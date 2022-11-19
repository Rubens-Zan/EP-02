#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

// GRR20206147 -> Rubens Zandomenighi Laszlo -> rzl20

int main ()
{
  // inicializa gerador de números aleatóreos
  srand(202202);
  
  SistLinear_t *SL;
  real_t *sol; 
  double tempo;
  real_t residue;

  // int diagTam[10] = {10, 30, 50, 128, 256, 512, 1000, 2000, 3000};
  SL = alocaSisLin(5);

  // for(int i=0;i < 9;++i){
    iniSisLin(SL,diagDominante, 10);
    
    eliminacaoGauss(SL,sol,&tempo);
    prnVetor(sol, SL->n);
    gaussSeidel(SL, sol,ERRO,&tempo);
    prnVetor(sol, SL->n); 

    refinamento(SL,sol,ERRO,&tempo);
    prnVetor(sol, SL->n); 

    normaL2Residuo(SL,sol,&residue);
    prnVetor(sol, SL->n); 

    // normaL2Residuo(SL,sol,)
    // liberaSisLin(SL); 
  // }

  // código do programa aqui
  
}

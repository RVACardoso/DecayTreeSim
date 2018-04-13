#include <iostream>

#ifndef __ODEpoint__
#define __ODEpoint__

//definição da struct que decreve um ponto num espaço de n dimensões
struct ODEpoint{
	double *var;// arrays geralmente usado para armazenar os tempos e as quantidades dos elementos 
	int size;//tamanho do array var
};

#endif

#include <iostream>
#include <cmath>
#include "element.h"
using namespace std;

element::element(): hltime(0.), spdecpb(0.), Nini(0.) {;}

element::element(double time1, double prob, double num): hltime(time1), spdecpb(prob), Nini(num) {;}

element::~element(){;}

//métodos para definição de elementos do private
void element::Sethltime(double time1){
	hltime=time1;
}

void element::Setspdecpb(double prob){
	spdecpb=prob;
}

// métodos para retornar a constante de decaimento de determinado elemento
double element::decaycnst(){
	return ((log(2.))/hltime);
}

//método para retornar a maior probabilidade das opções de decaimento
double element::bigdecprob(){
	return spdecpb;
}

//método para retornar a menor probabilidade das opções de decaimento
double element::smldecprob(){
	return (1.- spdecpb);
}
//metodo para retornar o número inicial de particulas de um determinado elemento
double element::getNini(){
	return Nini;
}

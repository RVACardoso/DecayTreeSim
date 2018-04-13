#include <iostream>
#include <vector>
#include "TFormula.h"
#include "ODEpoint.h"
using namespace std;

#ifndef __ODEsolver__
#define __ODEsolver__

class ODEsolver{
	public:
    	ODEsolver();
    	ODEsolver(vector<TFormula> vec);
    	void SetODEfunc(vector<TFormula>);
    	~ODEsolver();
    	//4 métodos distintos para resolução de sistemas de equações diferenciais
    	vector<ODEpoint> Euler(ODEpoint ini, double step, double T);
    	vector<ODEpoint> PredictorCorrector(ODEpoint ini, double step, double T);
    	vector<ODEpoint> RK2(ODEpoint ini, double step, double T);
    	vector<ODEpoint> RK4(ODEpoint ini, double step, double T);
	private:
	    vector<TFormula> F;	// vector com o sistema de equações a resolver
};

#endif
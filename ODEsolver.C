#include "ODEsolver.h"
#include "TFormula.h"
#include <vector>
#include <iostream>
using namespace std;

ODEsolver::ODEsolver(){;}

ODEsolver::ODEsolver(vector<TFormula> vec): F(vec) {;}

//definição das equações a resolver 
void ODEsolver::SetODEfunc(vector<TFormula> vec){
	F=vec;
}

ODEsolver::~ODEsolver(){;}

//aplicação do método de Euler
vector<ODEpoint> ODEsolver::Euler(ODEpoint ini, double step, double T){
	vector<ODEpoint > sol;
 	sol.push_back(ini);

	int funnb=F.size()+1;
	for(int j=0; sol[j].var[0] < T; ++j){
		sol.push_back(ODEpoint());
		sol[j+1].var = new double[ini.size];

		for (int i = 1; i < funnb; ++i){
			//calculo da função do ponto e obtenção do ponto seguinte
			sol[j+1].var[i] = sol[j].var[i] + step*F[i-1].EvalPar(sol[j].var);
		}
		//incremento do tempo
		sol[j+1].var[0]= sol[j].var[0]+step;
	}	
	return sol;
}

vector<ODEpoint> ODEsolver::PredictorCorrector(ODEpoint ini, double step, double T){
	vector<ODEpoint > sol;
 	sol.push_back(ini);

	for(int j=0; sol[j].var[0] < T; ++j){
		sol.push_back(ODEpoint());
		sol[j+1].var = new double[ini.size];	

		int funnb=F.size()+1;
		double F0, F1; 

		for (int i = 1; i < funnb; ++i){
			//calcular valor da função no ponto
			F0 = F[i-1].EvalPar(sol[j].var);

			//calcular valor da função no ponto seguinte
			double *temp=new double[ini.size];
			temp[0] = sol[j].var[0]+step;
			for (int k = 1; k < ini.size; ++k){
				temp[k]=sol[j].var[k] + step*F0;
			}
			F1 = F[i-1].EvalPar(temp);
			//determinção do ponto seguinte
			sol[j+1].var[i] = sol[j].var[i] + (step/2.)*(F0+F1);
		}
		//incremento do tempo
		sol[j+1].var[0]= sol[j].var[0]+step;
		//delete[] temp;
	}	
	//delete[] temp;
	return sol;
}
//aplicação do método Runge-Kutta 2
vector<ODEpoint> ODEsolver::RK2(ODEpoint ini, double step, double T){
	vector<ODEpoint > sol;
 	sol.push_back(ini);

	vector<double > k1;
	vector<double > k2;

	int funnb=F.size();
	int varnb=ini.size;
	for(int j=0; sol[j].var[0] < T; ++j){
		sol.push_back(ODEpoint());
		sol[j+1].var = new double[varnb];		

		//determinar K1's
		for (int i = 0; i < funnb; ++i){
			k1.push_back(double());
			k1[i]=step*F[i].EvalPar(sol[j].var); 
		}
		//determinar K2's
		double *temp;
		temp = new double[varnb];
		temp[0]=sol[j].var[0] + (step/2.);
		for (int i = 1; i < varnb; ++i){
			temp[i]= sol[j].var[i]+(1./2.)*k1[i-1];
		}
		for (int i = 0; i < funnb; ++i){
			k2.push_back(double());
			k2[i]=step*F[i].EvalPar(temp); 
		}
		
		//determinar proximos valor das variaveis
		for (int i = 1; i < funnb+1; ++i){
			sol[j+1].var[i] = sol[j].var[i] +k2[i-1];
		}
		//incremento do tempo
		sol[j+1].var[0]= sol[j].var[0]+step;
	}	
	return sol;
}

//aplicação do método Runge-Kutta 4
vector<ODEpoint> ODEsolver::RK4(ODEpoint ini, double step, double T){
	vector<ODEpoint > sol;
 	sol.push_back(ini);

	vector<double > k1,k2,k3,k4;

	int funnb=F.size();
	int varnb=ini.size;
	for(int j=0; sol[j].var[0] < T; ++j){
		sol.push_back(ODEpoint());
		sol[j+1].var = new double[ini.size];	

		//determinar K1's
		for (int i = 0; i < funnb; ++i){
			k1.push_back(double());
			k1[i]=F[i].EvalPar(sol[j].var); 
		}
		//determinar K2's
		double *temp = new double[varnb];
		temp[0]=sol[j].var[0]+(step/2.);
		for (int i = 1; i < varnb; ++i){
			temp[i]=sol[j].var[i] +(step/2.)*k1[i-1];
		}
		for (int i = 0; i < funnb; ++i){
			k2.push_back(double());
			k2[i]=F[i].EvalPar(temp); 
		}
		//determinar K3's 
		for (int i = 1; i < varnb; ++i){
			temp[i]=sol[j].var[i] + (step/2.)*k2[i-1];
		}
		for (int i = 0; i < funnb; ++i){
			k3.push_back(double());
			k3[i]=F[i].EvalPar(temp); 
		}
		//determinmar k4's
		temp[0]=sol[j].var[0]+step;

		for (int i = 1; i < varnb; ++i){
			temp[i]= sol[j].var[i] + step*k3[i-1];
		}
		for (int i = 0; i < funnb; ++i){
			k4.push_back(double());
			k4[i]=F[i].EvalPar(temp); 
		}
		//determinar proximos valor das variaveis
		for (int i = 1; i < funnb+1; ++i){
			sol[j+1].var[i] = sol[j].var[i] + (step/6.)*(k1[i-1]+2*k2[i-1]+2*k3[i-1]+k4[i-1]);
		}
		sol[j+1].var[0] = sol[j].var[0] + step;
	}	

	return sol;
}

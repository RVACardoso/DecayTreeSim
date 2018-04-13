#include "ODEdecay.h"
#include "MCdecay.h"
#include <string>
#include <iostream>
#include "ODEpoint.h"
#include <vector>
using namespace std;

//prototyping de funções definidas em Decayfunc.C
void Draw(double dt1, double dt2, double dt3, vector<ODEpoint > ode, vector<ODEpoint > mcbrt, vector<ODEpoint > mcpdf, int esc);
vector<double> Getmax(vector<ODEpoint > elem, int num);

int  main(){
	 // definir periodo de simulação
	double tmax=792.;//792h correspondem a 33 dias

	// definir propriedades e condições iniciais do elementos a serem simulados
	vector<element > elem;
	element Rn(14.61, 0.726 ,(2.85e15));
	elem.push_back(Rn);
	element At(7.21, 0.582, 0.);
	elem.push_back(At);
	element Po207(5.83, 1., 0.);
	elem.push_back(Po207);
	element Po211(1.43e-4, 1., 0.);
	elem.push_back(Po211);
	element Bi(288904.8, 1., 0.);
	elem.push_back(Bi);
	element Pb(1., 1., 0.); 
	elem.push_back(Pb);

	// ajuste do dt e simulação
	double dt1=0.0002;
	ODEdecay difdecay(dt1, tmax, elem);
	vector<ODEpoint > sol1;
	sol1=difdecay.decay();

	double dt2=0.1;
  	MCdecay mcdecaybrt(dt2, tmax, elem);
	vector<ODEpoint > sol2;
	sol2=mcdecaybrt.DecayBrt();

	// obtenção do máximo do 211Po através da função Getmax
	vector<double> max;
	max= Getmax(sol1, 5);
	cout << endl;
	cout << endl;
	cout << "Máximo do 211Po:  t-> " << max[0]<< "  nº-> "<< max[1]<< endl;
	cout << endl;

	//obtenção da quantidade final do 207Po 
	cout << "Passados 33 dias, a quantidade de 207Po é : " << endl;
	cout << "De acordo com a resolução da equação diferencial: " << endl;
	cout << " t: " <<sol1[(int)(tmax/dt1)].var[0] << " nº: " << sol1[(int)(tmax/dt1)].var[2] << endl;
	cout << "De acordo com técnicas de Monte-Carlo: " << endl;
	cout << " t: " << sol2[(int)(tmax/dt2-1.)].var[0] << " nº: " << sol2[(int)(tmax/dt2)].var[2] << endl;	
	cout << endl;

	//desenho dos gráficos que contém a evolução de todos os elementos para cada  um dos métodos
	difdecay.Draw();
	mcdecaybrt.Draw(0);
	// desenho dos gráficos de compração entre métodos, para cada um dos elementos
	Draw(dt1, dt2, 1., sol1, sol2, sol1, 0);

	// opção de simulação usando distrib. de probabilidade
	string res;
	while(1){
		cout << " Deseja simular o sistema usando técnicas de Monte-Carlo baseadas em distribuições de probabilidade, por um periodo de 100h? [y/n]" << endl;
		cin >> res;
		if(res=="y"){
			//realizar a simulação e desenhar resultados
			tmax=100.;
			double dt3=0.005;
			MCdecay mcdecaypdf(dt3, tmax, elem);
			vector<ODEpoint > sol3;
			sol3=mcdecaypdf.DecayPdf();
			mcdecaypdf.Draw(1);
			Draw(dt1, dt2, dt3, sol1, sol2, sol3, 1);
			break;
		}
		else if(res=="n"){
			break;
		}
		else
			cout << "Por favor, insira y ou n." << endl;
	}
	cout << endl;
}

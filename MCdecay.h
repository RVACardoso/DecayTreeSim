#ifndef __MCdecay__
#define __MCdecay__

#include <string>
#include "element.h"
#include "TGraph.h"
#include "ODEpoint.h"
#include <TFormula.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>    
#include <time.h>
using namespace std;

class MCdecay {
  public:
	MCdecay();
	MCdecay(double, double, vector<element > );
	~MCdecay();
	void setTmax(double);
	void setdt(double);
	vector<ODEpoint > DecayBrt();
	vector<ODEpoint > DecayPdf();
	void Draw(int n);
  private:
	double *tmp; //array com os instantes simulados
	double Tmax, dt;
	double *Rn, *At, *Bi, *Po207, *Po211, *Pb;//arrays com a quantidade dos elementos
	vector<element > elem;//vetor com as condiçoes inicias dos diferentes elementos
	vector<ODEpoint > sol;// resultado final;
};

#endif


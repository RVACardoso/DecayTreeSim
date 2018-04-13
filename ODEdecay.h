#include "TFormula.h"
#include "element.h"
#include "ODEsolver.h"
#include <iostream>
#include <vector>
#include "TGraph.h"

class ODEdecay{
	public:
		ODEdecay();
		ODEdecay(double intt, double ttot, vector<element > vec);
		vector<ODEpoint >  decay();
		void Draw();
		~ODEdecay();
	private:
		double dt;
		double tmax;
		vector<element > elem;//vetor com as condiçoes inicias dos diferentes elementos
		vector<ODEpoint > sol;// resultado final;
		
};
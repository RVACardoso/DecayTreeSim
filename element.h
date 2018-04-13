#include <iostream>
using namespace std;

#ifndef __element__
#define __element__

class element {
	public:
		element();
		element(double hltime, double spdecpb, double num);
		void Sethltime(double time1);
		void Setspdecpb(double prob);
		double getNini();
		~element();

		double decaycnst();
		double bigdecprob();
		double smldecprob();
	private:
		double hltime;//tempo de semi-vida do elemento
		double Nini;//numero inicial de particulas do elemento
		double spdecpb; //maior probabilidade das opçoes de decaimento do elemento
};

#endif
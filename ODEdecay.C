#include "ODEdecay.h"
#include <TAxis.h>
#include <TCanvas.h>
#include "TROOT.h"
#include <TFrame.h>
#include <cmath>
#include "TApplication.h"

#include <TImage.h>

using namespace std;

ODEdecay::ODEdecay():dt(0.), tmax(0.){;}

ODEdecay::ODEdecay(double intt, double ttot, vector<element > vec): dt(intt), tmax(ttot), elem(vec){;}

ODEdecay::~ODEdecay(){;}

//método para simulção do sistema recorrendo à resolução numerica de um sistema de equações diferenciais
vector<ODEpoint > ODEdecay::decay(){

	//criação das equaçoes diferenciais que descrevem o comportamento de cada um dos elementos
 	TFormula  s1("s1",    "-[0]*x[1]");//rn
    s1.SetParameter(0, elem[0].decaycnst());  
  	TFormula  s2("s2", "-[2]*x[2]+[1]*[0]*x[1]"); //Po207
  	s2.SetParameter(0, elem[0].decaycnst());  
    s2.SetParameter(1, elem[0].smldecprob()); 
    s2.SetParameter(2, elem[2].decaycnst()); 
    TFormula  s3("s3", "-[3]*x[3]+[0]*x[2]+[1]*[2]*x[4]");//Bi
    s3.SetParameter(0, elem[2].decaycnst()); 
    s3.SetParameter(1, elem[1].smldecprob()); 
    s3.SetParameter(2, elem[1].decaycnst()); 
    s3.SetParameter(3, elem[4].decaycnst());  
    TFormula  s4("s3", "-[2]*x[4]+[0]*[1]*x[1]");//At
    s4.SetParameter(0, elem[0].bigdecprob());
    s4.SetParameter(1, elem[0].decaycnst());
    s4.SetParameter(2, elem[1].decaycnst()); 
    TFormula  s5("s5", "-[2]*x[5]+[0]*[1]*x[4]");//Po211
    s5.SetParameter(0, elem[1].bigdecprob());  
    s5.SetParameter(1, elem[1].decaycnst());
    s5.SetParameter(2, elem[3].decaycnst()); 
    TFormula  s6("s6","[1]*x[3]+[0]*x[5]"); //Pb
    s6.SetParameter(0, elem[3].decaycnst());
    s6.SetParameter(1, elem[4].decaycnst());

  	//criação de vetor com as equações
  	vector<TFormula> vec;
  	vec.push_back(s1);
  	vec.push_back(s2);
  	vec.push_back(s3);
  	vec.push_back(s4);
  	vec.push_back(s5);
  	vec.push_back(s6);
    ODEsolver solv(vec);

    //condiçoes inicias
    ODEpoint ini;
    ini.size=7.;
    ini.var = new double[ini.size];
    ini.var[0]=0.;//tempo inicial
    ini.var[1] = elem[0].getNini();
    ini.var[2] = elem[1].getNini();
    ini.var[3] = elem[2].getNini();
    ini.var[4] = elem[3].getNini();
    ini.var[5] = elem[4].getNini();
    ini.var[6] = elem[5].getNini();

    //resolução do sistema utilizando Runge-Kutta 4
    sol=solv.RK4(ini,dt,tmax); 
    for(int i=0; i<sol.size();++i){
    	sol[i].var[1] = round(sol[i].var[1]);
    	sol[i].var[2] = round(sol[i].var[2]);
    	sol[i].var[3] = round(sol[i].var[3]);
    	sol[i].var[4] = round(sol[i].var[4]);
    	sol[i].var[5] = round(sol[i].var[5]);
    	sol[i].var[6] = round(sol[i].var[6]);
	}
    return sol;
}

//método para desenhar o gráfico com as quantidades de todos os elementos obtidas pela resolução das equaçoes
void ODEdecay::Draw(){
// alocação de memórias
 	int arrsize=(100./dt +1.);
    double *tmp = new double[arrsize];
    double *Rn = new double[arrsize];
    double *At = new double[arrsize];
    double *Po211 = new double[arrsize];
    double *Po207 = new double[arrsize];
    double *Bi = new double[arrsize];
    double *Pb = new double[arrsize];


    for(int i=0; i< arrsize; ++i){
	    tmp[i] = sol[i].var[0];
	    Rn[i] = sol[i].var[1];
	    At[i] = sol[i].var[4];
	    Po207[i] = sol[i].var[2];
	    Bi[i] = sol[i].var[3];
	    Po211[i] = sol[i].var[5];
	    Pb[i] = sol[i].var[6];
    }

	// criação dos objetos necessários à utilização do root
	TApplication *gMyRootApp;
 	int argc = 0;
  	char **argv = NULL;
 	gMyRootApp = new TApplication("Comput Phys Application", &argc, argv);
  	gMyRootApp->SetReturnFromRun(true);
  	gROOT->SetStyle("Plain");

	//criação e desenho dos canvas e dos gráficos necessários
	TCanvas *c1 = new TCanvas("c1","Evolução do sistema");//,200,10,700,500);

   	c1->SetFillColor(42);
   	c1->SetGrid();

   	TGraph *g1 = new TGraph(arrsize, tmp , Rn); //Rn
 	g1->SetMarkerStyle(20);
  	g1->SetMarkerColor(kRed);
	g1->SetMarkerSize(0.5);
  	g1->SetTitle("ODE method");
  	g1->GetXaxis()->SetTitle("Time(h)");
    g1->GetYaxis()->SetTitle("Number of particles");
    g1->GetYaxis()->SetTitleOffset(1.2);
    g1->GetXaxis()->CenterTitle();
    g1->GetYaxis()->CenterTitle();
    g1->Draw("AP"); 
    TGraph *g2 = new TGraph(arrsize,tmp,At); //At
    g2->SetMarkerStyle(20);
    g2->SetMarkerColor(kBlue);
	g2->SetMarkerSize(0.5);
	g2->Draw("P"); 
	TGraph *g3 = new TGraph(arrsize,tmp,Po207); //Po207
    g3->SetMarkerStyle(20);
    g3->SetMarkerColor(kGreen);
	g3->SetMarkerSize(0.5);
	g3->Draw("P"); 
	TGraph *g4 = new TGraph(arrsize,tmp,Po211); //Po211
    g4->SetMarkerStyle(20);
    g4->SetMarkerColor(kYellow);
	g4->SetMarkerSize(0.5);
	g4->Draw("P"); 
	TGraph *g5= new TGraph(arrsize,tmp,Bi); //Bi 
    g5->SetMarkerStyle(20);
    g5->SetMarkerColor(kViolet);
	g5->SetMarkerSize(0.5);
	g5->Draw("P"); 
	TGraph *g6 = new TGraph(arrsize,tmp, Pb); //Pb
    g6->SetMarkerStyle(20);
    g6->SetMarkerColor(kBlack);
	g6->SetMarkerSize(0.5);
	g6->Draw("P"); 
 
    c1->Update();
    c1->Modified();
    c1->WaitPrimitive();
 	TImage *img = TImage::Create();
   	img->FromPad(c1);
	img->WriteImage("ODE.png");

  

   //desalocação de memória
	delete[] tmp; 
	delete[] Rn; 
	delete[] At; 
	delete[] Po211;
	delete[] Po207;
	delete[] Bi; 
	delete[] Pb;
}

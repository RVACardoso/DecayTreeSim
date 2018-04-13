#include "MCdecay.h" 
#include <TRandom3.h>
#include <TAxis.h>
#include <TCanvas.h>
#include "TROOT.h"
#include <TFrame.h>
#include <cmath>
#include "TApplication.h"
#include <TImage.h>

TRandom3 Rnd(time(NULL));

MCdecay::MCdecay(){
	Tmax=0;
	dt=0;
	tmp=NULL;
	Rn = NULL;
	At = NULL;
	Bi = NULL;
	Po211 = NULL;
	Po207 = NULL;
	Pb = NULL;
}

MCdecay::MCdecay(double intt, double tmax, vector<element > vec ): Tmax(tmax), dt(intt), elem(vec){
	int arrsize=(Tmax/dt)+1;
	tmp= new double[arrsize];
	Rn = new double[arrsize];
	At = new double[arrsize];
	Bi = new double[arrsize];
	Po211 = new double[arrsize];
	Po207 = new double[arrsize];
	Pb = new double[arrsize];
}

MCdecay::~MCdecay(){
	delete[] tmp;
	delete[] Rn;
	delete[] At;
	delete[] Bi;
	delete[] Po211;
	delete[] Po207;
	delete[] Pb;
}

//método para definição de elementos do private
void MCdecay::setTmax(double tmax){
	Tmax=tmax;
}

void MCdecay::setdt(double intt){
	dt=intt;
}

//simulação usando técnica de Monte-Carlo na qual se simulam ~10e5 particulas e se extrapola o resultado às restantes
vector<ODEpoint > MCdecay::DecayBrt(){
	//inicialização das variaveis
	double T=0., R;
	Rn[0]=elem[0].getNini(), At[0]=elem[1].getNini(), Bi[0]=elem[4].getNini();
	Po207[0]=elem[2].getNini(), Po211[0]=elem[3].getNini(), Pb[0]=elem[5].getNini(); 	

	double NRn=Rn[0], NAt=At[0], NBi=Bi[0], NPo207=Po207[0], NPo211=Po211[0], NPb=Pb[0];

	double redfac=1e10;	// fator de redução que seleciona a porção de particulas a simular
	NRn=NRn/redfac;

	int arrsize=(Tmax/dt)+1;
	int i=0;
	//loop ao longo do tempo para avaliar evolução temporal do sistema
	while (T < Tmax){
		if(i%(arrsize/10)==0)
			cout << " #"<< flush;
		//loops ao longo de cada uma das particulas de determinado elemento
		for (int n = 0; n < NRn; ++n){
			R =  Rnd.Uniform(0,1);//geração de aleatorio cujo valor determina se a particula decai ou nao
			if (R< elem[0].decaycnst()*dt){// se o valor for inferior ao produto do dt pela constante de decaimento, há decaimento
				--NRn; 
				if(R<=0.274*elem[0].decaycnst()*dt){//determinação, com base no aleatório já gerado, da opção de decaimento
					++NPo207;
				}
				else{
					++NAt;
				}
			}
		}
		for (int n = 0; n < NAt; ++n){
			R =  Rnd.Uniform(0,1);
			if (R< elem[1].decaycnst()*dt){
				--NAt;
				if(R<=0.418*elem[1].decaycnst()*dt)
					++NBi;
				else
					++NPo211;
			}
		}
		for (int n = 0; n < NPo207; ++n){
			R =  Rnd.Uniform(0,1);
			if (R< elem[2].decaycnst()*dt){
				--NPo207;
				++NBi;
			}
		}
		for (int n = 0; n < NBi; ++n){
			R =  Rnd.Uniform(0,1);
			if (R< elem[4].decaycnst()*dt){
				--NBi;
				++NPb;
			}
		}
		for (int n = 0; n <= NPo211; ++n){
			R =  Rnd.Uniform(0,1);
			if (R< elem[3].decaycnst()*dt){
				--NPo211;
				++NPb;
			}
		}

	//armazenamento dos valores atuais depois de haver extrapolação para a totalidade das particulas
		T+=dt;
		tmp[i]=T;
		Rn[i]=NRn*redfac;
		At[i]=NAt*redfac;
		Bi[i]=NBi*redfac;
		Po207[i]=NPo207*redfac;
		Po211[i]=NPo211*redfac/(500);
		Pb[i]=NPb*redfac;
		++i;
	}
// criação, recorrendo aos dados disponiveis, de um elemento de ODEpoint que atua como resultado da simulação
	for (int i = 0; i < arrsize; ++i){
		ODEpoint psol;
		psol.size=7;
		psol.var = new double[psol.size];
		psol.var[0]=tmp[i];
		psol.var[1]=Rn[i];
		psol.var[2]=At[i];
		psol.var[3]=Po207[i];
		psol.var[4]=Po211[i];
		psol.var[5]=Bi[i];
		psol.var[6]=Pb[i];
		sol.push_back(psol);
	}
	return sol;
}

//funçaõ de acceptance-rejection para geração de numeros aleatorios entre 0 e 1 de acordo com uma distrib. de probabilidade normal
//esta função é chamada abaixo 
double randgen(double y,double o){// y é a média, o é a variância
	if(y==0)
		return 0;
	else{
		TFormula *pdf = new TFormula("pdf","(1/z*sqrt(2*TMath::Pi()))*exp(-((x-y)^2/(2*z*z)))"); 
		double pmax= pdf->Eval(y,y,o);
		while(1){
			double R1 = Rnd.Uniform(0,1);
			double pxr = pdf->Eval(R1,y,o);
			double R2 = Rnd.Uniform(0,1);
			if(R2 <= pxr/pmax){
				return R1;
			}
		}
	}
}

//simulação usando técnicas de Monte-Carlo com base em distribuições de probabilidade
vector<ODEpoint > MCdecay::DecayPdf(){
	//inicialização das variáveis
	double T=0., R,R2;

	Rn[0]=elem[0].getNini(), At[0]=elem[1].getNini(), Bi[0]=elem[4].getNini();
	Po207[0]=elem[2].getNini(), Po211[0]=elem[3].getNini(), Pb[0]=elem[5].getNini(); 	
	double NRn=Rn[0], NAt=At[0], NBi=Bi[0], NPo207=Po207[0], NPo211=Po211[0], NPb=Pb[0];

	int arrsize=(Tmax/dt)+1;
	double varRn, varAt, varPo207, varPo211, varBi;

	int i=0;
	double o=0.5;
	//loop ao longo do tempo para simulação da evolução temporal
	while (T < Tmax){
		if(i%(arrsize/10)==0)
			cout << " #"<< flush;
		//geração, usando a função randgen, da fração de particulas de determinado elemento que decaem
		//determinação da variação do numero de particulas com base nesta fração
		varRn=round(randgen(elem[0].decaycnst()*dt, elem[0].decaycnst()*dt*o)*NRn);
		NRn-=varRn;
		if(Rnd.Uniform(0,1) <=0.274)// determinação da opção de decaimento
			NPo207+=varRn;
		else
			NAt+=varRn;

		varAt=round(randgen(elem[1].decaycnst()*dt,elem[1].decaycnst()*dt*o)*NAt);
		NAt-=varAt;
		if(Rnd.Uniform(0,1) <=0.418)
			NBi+=varAt;
		else
			NPo211+=varAt;

		varPo207=round(randgen(elem[2].decaycnst()*dt,elem[2].decaycnst()*dt*o)*NPo207);
		NPo207-=varPo207;
		NBi+=varPo207;

		/*
		De forma a aumentar bastante a eficácia computacional, ignorou-se o decaimento do Bismuto já que o 
		seu periodo de semi-vida é tão elevado

		varBi=round(randgen(elem[4].decaycnst()*dt)*NBi); //elem[4]
		NBi-=varBi;
		NPb+=varBi;*/

		varPo211=round(randgen(elem[3].decaycnst()*dt, elem[3].decaycnst()*dt*o)*NPo211);
		NPo211-=varPo211;
		NPb+=varPo211;

		//armazenamento das quantidade atuais dos elementos
		T+=dt;
		tmp[i]=T;
		Rn[i]=NRn;
		At[i]=NAt;
		Bi[i]=NBi;
		Po207[i]=NPo207;
		Po211[i]=NPo211;
		Pb[i]=NPb;
		++i;
		
	}

	// criação, recorrendo aos dados disponiveis, de um elemento de ODEpoint que atua como resultado da simulação
	for (int i = 0; i < arrsize; ++i){
		ODEpoint psol;
		psol.size=7;
		psol.var = new double[psol.size];
		psol.var[0]=tmp[i];
		psol.var[1]=Rn[i];
		psol.var[2]=At[i];
		psol.var[3]=Po207[i];
		psol.var[4]=Po211[i];
		psol.var[5]=Bi[i];
		psol.var[6]=Pb[i];
		sol.push_back(psol);
	}
	return sol;
}

//método para desenhar os graficos dos elementos simulados por tecnicas de Monte-Carlo
void MCdecay::Draw(int n){
	//determinção dos tamanho dos arrays e alocação de memória
	int arrsize=(100./dt)+1;

	//criação e desenho dos canvas e graphs necessários
	TCanvas *c2 = new TCanvas("c2","Evolução do sistema");//,200,10,700,500);
   	c2->SetFillColor(42);
   	c2->SetGrid();
		
 
  	TGraph *g1 = new TGraph(arrsize, tmp,Rn); 
  	g1->SetMarkerStyle(20);
  	g1->SetMarkerColor(kRed);
  	g1->SetMarkerSize(0.5);
  	if(n==0)
  		g1->SetTitle("Monte-Carlo Brt method");
  	else	
  		g1->SetTitle("Monte-Carlo PDF method");
  	g1->GetXaxis()->SetTitle("Time(h)");
  	g1->GetYaxis()->SetTitle("Number of particles");
  	g1->GetYaxis()->SetTitleOffset(1.2);
  	g1->GetXaxis()->CenterTitle();
  	g1->GetYaxis()->CenterTitle();
  	g1->Draw("AP"); 
  	TGraph *g2 = new TGraph(arrsize,tmp,At); 
    g2->SetMarkerStyle(20);
    g2->SetMarkerColor(kBlue);
    g2->SetMarkerSize(0.5);
	g2->Draw("P"); 
	TGraph *g3 = new TGraph(arrsize,tmp,Po207); 
	g3->SetMarkerStyle(20);
	g3->SetMarkerColor(kGreen);
	g3->SetMarkerSize(0.5);
	g3->Draw("P"); 
	TGraph *g4 = new TGraph(arrsize,tmp,Po211); 
	g4->SetMarkerStyle(20);
	g4->SetMarkerColor(kYellow);
	g4->SetMarkerSize(0.5);
	g4->Draw("P"); 
	TGraph *g5= new TGraph(arrsize,tmp,Bi); 
	g5->SetMarkerStyle(20);
	g5->SetMarkerColor(kViolet);
	g5->SetMarkerSize(0.5);
	g5->Draw("P"); 
	TGraph *g6 = new TGraph(arrsize,tmp,Pb); 
	g6->SetMarkerStyle(20);
	g6->SetMarkerColor(kBlack);
	g6->SetMarkerSize(0.5);
	g6->Draw("P"); 

	c2->Update();
	c2->Modified();
 	c2->WaitPrimitive();
 	
	TImage *img = TImage::Create();
   	img->FromPad(c2);
	
	if(n==0)
  		img->WriteImage("MCbrt.png");
  	else	
  		img->WriteImage("MCpdf.png");
 

  	
 	
  	//desalocação da memória usada
 	delete c2;
}




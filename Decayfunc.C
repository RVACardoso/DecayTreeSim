#include "TGraph.h"
#include <iostream>
#include "ODEpoint.h"
#include <TAxis.h>
#include <TCanvas.h>
#include "TROOT.h"
#include <TFrame.h>
#include <cmath>
#include "TApplication.h"
#include <TImage.h>

using namespace std;

//função de desenho dos gráfico que mostram, para cada elemento, a comparação entre os dois/três métodos
void Draw(double dt1, double dt2, double dt3, vector<ODEpoint > ode, vector<ODEpoint > mcbrt, vector<ODEpoint > mcpdf, int num){
  //determinar tamanho dos arrays a utilizar
  int arrsize1=(100./dt1)+1.;
  int arrsize2=(100./dt2)+1.;
  int arrsize3=(100./dt3)+1.;

  //alocar dinamicamente o espaço necessário para os arrays
  double *tmp1 =new double[arrsize1];
  double *Rn1 =new double[arrsize1];
  double *At1 =new double[arrsize1];
  double *Po2071=new double[arrsize1];
  double *Po2111=new double[arrsize1];
  double *Bi1=new double[arrsize1];
  double *Pb1=new double[arrsize1];

  double *tmp2 =new double[arrsize2];
  double *Rn2 =new double[arrsize2];
  double *At2 =new double[arrsize2];
  double *Po2072=new double[arrsize2];
  double *Po2112=new double[arrsize2];
  double *Bi2=new double[arrsize2];
  double *Pb2=new double[arrsize2];

  double *tmp3 =new double[arrsize3];
  double *Rn3 =new double[arrsize3];
  double *At3 =new double[arrsize3];
  double *Po2073=new double[arrsize3];
  double *Po2113=new double[arrsize3];
  double *Bi3=new double[arrsize3];
  double *Pb3=new double[arrsize3];

	//converter os vectors de ODEpoint em arrays apropriados ao desenho do gráfico
  for (int i = 0; i < arrsize1; ++i){
    tmp1[i]=ode[i].var[0];
    Rn1[i]=ode[i].var[1];
    At1[i]=ode[i].var[4];
    Po2071[i]=ode[i].var[2];
    Po2111[i]=ode[i].var[5];
    Bi1[i]=ode[i].var[3];
    Pb1[i]=ode[i].var[6];
    }
  for (int i = 0; i < arrsize2; ++i){
    tmp2[i]=mcbrt[i].var[0];
    Rn2[i]=mcbrt[i].var[1];
    At2[i]=mcbrt[i].var[2];
    Po2072[i]=mcbrt[i].var[3];
    Po2112[i]=mcbrt[i].var[4];
    Bi2[i]=mcbrt[i].var[5];
    Pb2[i]=mcbrt[i].var[6];
  }
  if(num==1){
    for (int i = 0; i < arrsize3; ++i){
      tmp3[i]=mcpdf[i].var[0];
      Rn3[i]=mcpdf[i].var[1];
      At3[i]=mcpdf[i].var[2];
      Po2073[i]=mcpdf[i].var[3];
      Po2113[i]=mcpdf[i].var[4];
      Bi3[i]=mcpdf[i].var[5];
      Pb3[i]=mcpdf[i].var[6];
    }
  }
//criação e desenho do canvas e dos graficos necessários
  TCanvas *c11 = new TCanvas("c11","Evolução do sistema");
  c11->SetFillColor(42);
  c11->SetGrid();
   
  TGraph *g11 = new TGraph(arrsize1,tmp1 ,Rn1); 
  g11->SetMarkerStyle(20);
  g11->SetMarkerColor(kRed);
  g11->SetMarkerSize(0.5);
  g11->SetTitle("Radon(211Rn)");
  g11->GetXaxis()->SetTitle("Time(h)");
  g11->GetYaxis()->SetTitle("Number of particles");
  g11->GetYaxis()->SetTitleOffset(1.2);
  g11->GetXaxis()->CenterTitle();
  g11->GetYaxis()->CenterTitle();
  g11->Draw("AP");
  TGraph *g12 = new TGraph(arrsize2,tmp2,Rn2); 
  g12->SetMarkerStyle(20);
  g12->SetMarkerColor(kBlue);
  g12->SetMarkerSize(0.5);
  g12->Draw("P");
  if(num==1){//condição para desenhar o grafico gerado usando distrib. de probabilidade apenas quando pedido
  	TGraph *g13 = new TGraph(arrsize3,tmp3,Rn3); 
  	g13->SetMarkerStyle(20);
 	g13->SetMarkerColor(kGreen);
 	g13->SetMarkerSize(0.5);
  	g13->Draw("P");
  }
  c11->Update();
  c11->Modified();
  c11->WaitPrimitive();
  TImage *img1 = TImage::Create();
  img1->FromPad(c11);
  img1->WriteImage("Radon.png");


  TCanvas *c12 = new TCanvas("c12","Evolução do sistema");
  c12->SetFillColor(42);
  c12->SetGrid();

  TGraph *g21 = new TGraph(arrsize1,tmp1 ,At1); 
  g21->SetMarkerStyle(20);
  g21->SetMarkerColor(kRed);
  g21->SetMarkerSize(0.5);
  g21->SetTitle("Astatine(211At)");
  g21->GetXaxis()->SetTitle("Time(h)");
  g21->GetYaxis()->SetTitle("Number of particles");
  g21->GetYaxis()->SetTitleOffset(1.2);
  g21->GetXaxis()->CenterTitle();
  g21->GetYaxis()->CenterTitle();
  g21->Draw("AP");
  TGraph *g22 = new TGraph(arrsize2,tmp2,At2); 
  g22->SetMarkerStyle(20);
  g22->SetMarkerColor(kBlue);
  g22->SetMarkerSize(0.5);
  g22->Draw("P");
  if(num==1){
  	TGraph *g23 = new TGraph(arrsize3,tmp3,At3); 
  	g23->SetMarkerStyle(20);
  	g23->SetMarkerColor(kGreen);
  	g23->SetMarkerSize(0.5);
  	g23->Draw("P");
  }

  c12->Update();
  c12->Modified();
  c12->WaitPrimitive();
  TImage *img2 = TImage::Create();
  img2->FromPad(c12);
  img2->WriteImage("Astatine.png");


  TCanvas *c13 = new TCanvas("c13","Evolução do sistema");
  c13->SetFillColor(42);
  c13->SetGrid();

  TGraph *g31 = new TGraph(arrsize1,tmp1 ,Po2071); 
  g31->SetMarkerStyle(20);
  g31->SetMarkerColor(kRed);
  g31->SetMarkerSize(0.5);
  g31->SetTitle("Polonium(207Po)");
  g31->GetXaxis()->SetTitle("Time(h)");
  g31->GetYaxis()->SetTitle("Number of particles");
  g31->GetYaxis()->SetTitleOffset(1.2);
  g31->GetXaxis()->CenterTitle();
  g31->GetYaxis()->CenterTitle();
  g31->Draw("AP"); 
  TGraph *g32 = new TGraph(arrsize2,tmp2, Po2072); 
  g32->SetMarkerStyle(20);
  g32->SetMarkerColor(kBlue);
  g32->SetMarkerSize(0.5);
  g32->Draw("P"); 
  if(num==1){
  	TGraph *g33 = new TGraph(arrsize3,tmp3,Po2073); 
  	g33->SetMarkerStyle(20);
  	g33->SetMarkerColor(kGreen);
  	g33->SetMarkerSize(0.5);
  	g33->Draw("P"); 
  }
  c13->Update();
  c13->Modified();
  c13->WaitPrimitive();
  TImage *img3 = TImage::Create();
  img3->FromPad(c13);
  img3->WriteImage("Polonium207.png");



  TCanvas *c14 = new TCanvas("c14","Evolução do sistema");//,200,10,700,500);
  c14->SetFillColor(42);
  c14->SetGrid();

  TGraph *g41 = new TGraph(arrsize1,tmp1 ,Bi1); 
  g41->SetMarkerStyle(20);
  g41->SetMarkerColor(kRed);
  g41->SetMarkerSize(0.5);
  g41->SetTitle("Bismuth(207Bi)");
  g41->GetXaxis()->SetTitle("Time(h)");
  g41->GetYaxis()->SetTitle("Number of particles");
  g41->GetYaxis()->SetTitleOffset(1.2);
  g41->GetXaxis()->CenterTitle();
  g41->GetYaxis()->CenterTitle();
  g41->Draw("AP"); 
  TGraph *g42 = new TGraph(arrsize2,tmp2, Bi2); 
  g42->SetMarkerStyle(20);
  g42->SetMarkerColor(kBlue);
  g42->SetMarkerSize(0.5);
  g42->Draw("P"); 
  if(num==1){
  	TGraph *g43 = new TGraph(arrsize3,tmp3,Bi3); 
  	g43->SetMarkerStyle(20);
  	g43->SetMarkerColor(kGreen);
  	g43->SetMarkerSize(0.5);
  	g43->Draw("P"); 
  }
  c14->Update();
  c14->Modified();
  c14->WaitPrimitive();
  TImage *img4 = TImage::Create();
  img4->FromPad(c14);
  img4->WriteImage("Bismuth.png");

  TCanvas *c15 = new TCanvas("c15","Evolução do sistema");//,200,10,700,500);
  c15->SetFillColor(42);
  c15->SetGrid();

  TGraph *g51 = new TGraph(arrsize1,tmp1 ,Po2111); 
  g51->SetMarkerStyle(20);
  g51->SetMarkerColor(kRed);
  g51->SetMarkerSize(0.5);
  g51->SetTitle("Polonium(211Po)");
  g51->GetXaxis()->SetTitle("Time(h)");
  g51->GetYaxis()->SetTitle("Number of particles");
  g51->GetYaxis()->SetTitleOffset(1.2);
  g51->GetXaxis()->CenterTitle();
  g51->GetYaxis()->CenterTitle();
  g51->Draw("AP"); 
  TGraph *g52 = new TGraph(arrsize2,tmp2, Po2112); 
  g52->SetMarkerStyle(20);
  g52->SetMarkerColor(kBlue);
  g52->SetMarkerSize(0.5);
  g52->Draw("P"); 
  if(num==1){
  	TGraph *g53 = new TGraph(arrsize3,tmp3,Po2113); 
  	g53->SetMarkerStyle(20);
  	g53->SetMarkerColor(kGreen);
  	g53->SetMarkerSize(0.5);
  	g53->Draw("P"); 
  }

  c15->Update();
  c15->Modified();
  c15->WaitPrimitive();
  TImage *img5 = TImage::Create();
  img5->FromPad(c15);
  img5->WriteImage("Polonium211.png");


  TCanvas *c16 = new TCanvas("c16","Evolução do sistema");//,200,10,700,500);
  c16->SetFillColor(42);
  c16->SetGrid();

  TGraph *g61 = new TGraph(arrsize1,tmp1 ,Pb1); 
  g61->SetMarkerStyle(20);
  g61->SetMarkerColor(kRed);
  g61->SetMarkerSize(0.5);
  g61->SetTitle("Lead(207Pb)");
  g61->GetXaxis()->SetTitle("Time(h)");
  g61->GetYaxis()->SetTitle("Number of particles");
  g61->GetYaxis()->SetTitleOffset(1.2);
  g61->GetXaxis()->CenterTitle();
  g61->GetYaxis()->CenterTitle();
  g61->Draw("AP"); 
  TGraph *g62 = new TGraph(arrsize2,tmp2, Pb2); 
  g62->SetMarkerStyle(20);
  g62->SetMarkerColor(kBlue);
  g62->SetMarkerSize(0.5);
  g62->Draw("P"); 
  if(num==1){
  	TGraph *g63 = new TGraph(arrsize3,tmp3,Pb3); 
  	g63->SetMarkerStyle(20);
  	g63->SetMarkerColor(kGreen);
  	g63->SetMarkerSize(0.5);
  	g63->Draw("P"); 
  }

  c16->Update();
  c16->Modified();
  c16->WaitPrimitive();
  TImage *img6 = TImage::Create();
  img6->FromPad(c16);
  img6->WriteImage("Lead.png");


	//desalocação da memória usada
  delete[] tmp1;
  delete[] Rn1; 
  delete[] At1; 
  delete[] Po2071;
  delete[] Po2111;
  delete[] Bi1;
  delete[] Pb1;

  delete[] tmp2; 
  delete[] Rn2; 
  delete[] At2; 
  delete[] Po2072;
  delete[] Po2112;
  delete[] Bi2;
  delete[] Pb2;
  if(num==1){
  	delete[] tmp3; 
  	delete[] Rn3; 
  	delete[] At3;
  	delete[] Po2073;
  	delete[] Po2113;
  	delete[] Bi3;
  	delete[] Pb3;
  }
  delete c11;
  delete c12;
  delete c13;
  delete c14;
  delete c15;
  delete c16;
}

//função para obter a quantidade máxima de um determinado elemento e o tempo em que esta ocorre
vector<double> Getmax(vector<ODEpoint > elem, int num){
	vector<double> max(2);
	double M,m;
	//procura em todos os pontos uma quantidade maior que a guardada em M, se a encontrar guarda-a
	for(int i=0; i< elem.size() ;++i){
		if(elem[i].var[num]>M){
			M=elem[i].var[num];
			m=i;
		}
	}
	max[0]=elem[m].var[0];
	max[1]=M;
	return max;
}
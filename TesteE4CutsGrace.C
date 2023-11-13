/*
Inv Mass
0: LO
1: NLO

With M_{recoil} € [120,135] yes or no
*/
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "TMath.h"
#include <vector>
#else
class ExRootTreeReader;
class ExRootResult;
#endif

template<typename T>
void Function_pushback(const TClonesArray& inColl ,vector<T*>& outColl)//, Double_t ptMin=15, Double_t ptMax= 70)
{

  const TObject *object;

  for (Int_t i = 0; i < inColl.GetEntriesFast(); i++)
  {

   object = inColl.At(i);
   const T *t = static_cast<const T*>(object);

   //if(t->P4().Pt() < ptMin && t->P4().Pt() > ptMax) continue;
   //if(TMath::Abs(t->P4().Eta()) > etaMax) continue;

   outColl.push_back(t);

  }
}
template<typename T1>
void BTagging(const TClonesArray& inColl , vector<T1*>& outColl)
{
    const TObject *objectb;

    for (Int_t i = 0; i < inColl.GetEntriesFast(); i++)
    {
        objectb = inColl.At(i);
        const T1 *t1 = static_cast<const T1*>(objectb);

        if(((t1->BTag) == 1) || ((t1->BTag) == 3)) outColl.push_back(t1);
    }
    
}
/*----- Function to Study the Kinematics and Fill the Histograms      */
void AnalyseEvents(ExRootTreeReader *treeReader, TH1F *hist1, Int_t EnergyCM,Int_t cuts,Int_t cutsee,Double_t sc)
{
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  
  Long64_t numberOfEntries = treeReader->GetEntries();

  const Electron *elec1, *elec2;
  const Jet *bjet1, *bjet2;
  const Muon *muon1,  *muon2;
  const Photon *photon1, *photon2;
  const MissingET *met1, *met2;

  vector<const Electron*> *Velectrons = new vector<const Electron*>();;
  vector<const Jet*> *Vbjets = new vector<const Jet*>();;
  vector<const Muon*> *Vmuons = new vector<const Muon*>();;
  vector<const Photon*> *Vphotons = new vector<const Photon*>();; 
  //vector<const MissingET> *Vmet = new vector<const MissingET*>();;

  Double_t MassBB, MassEE,EnergyEE, RecoilMass, Met;
  Double_t dre1j1, dre1j2, dre2j1,dre2j2, drj1j2;
  Double_t drmu1j1, drmu1j2, drmu2j1, drmu2j2;
  Double_t drpho1j1, drpho1j2, drpho2j1, drpho2j2;
  Double_t Energy4;
  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
    Velectrons->clear();
    Vbjets->clear();
    Vmuons->clear();
    Vphotons->clear();

    Function_pushback(*branchElectron, *Velectrons);
    Function_pushback(*branchMuon, *Vmuons);
    Function_pushback(*branchPhoton, *Vphotons);
    BTagging(*branchJet, *Vbjets);
    

    if(Velectrons->size() < 2) continue;
    if(Vmuons->size() < 2) continue;
    if(Vphotons->size() < 2) continue;
    if(Vbjets->size() < 2) continue;

    elec1 = Velectrons->at(0);
    elec2 = Velectrons->at(1);

    muon1 = Vmuons->at(0);
    muon2 = Vmuons->at(1);

    bjet1 = Vbjets->at(0);
    bjet2 = Vbjets->at(1);
  
    photon1 = Vphotons->at(0);
    photon2 = Vphotons->at(1);
 /*
    elec1->Charge == elec2 ->Charge) continue;
   
    if(muon1->Charge == muon2 ->Charge) continue;
    
//  I wnat my jets to have a DeltaR min of 0.5 between electrons and jets 
    dre1j1 = (elec1->P4()).DeltaR(bjet1->P4());
    dre1j2 = (elec1->P4()).DeltaR(bjet2->P4());
    dre2j1 = (elec2->P4()).DeltaR(bjet1->P4());
    dre2j2 = (elec2->P4()).DeltaR(bjet2->P4());
    drj1j2 = (bjet1->P4()).DeltaR(bjet2->P4());
    drmu1j1 = (muon1->P4()).DeltaR(bjet1->P4());
    drmu1j2 = (muon1->P4()).DeltaR(bjet2->P4());
    drmu2j1 = (muon2->P4()).DeltaR(bjet1->P4());
    drmu2j2 = (muon2->P4()).DeltaR(bjet2->P4());
    drpho1j1 = (photon1->P4()).DeltaR(bjet1->P4());
    drpho1j2 = (photon1->P4()).DeltaR(bjet2->P4());
    drpho2j1 = (photon2->P4()).DeltaR(bjet1->P4());
    drpho2j2 = (photon2->P4()).DeltaR(bjet2->P4());

    if ((dre1j1 < 0.5) || (dre1j2 < 0.5) || (dre2j1 < 0.5) || (dre2j2 < 0.5) || (drmu1j1 < 0.5) || (drmu1j2 <  0.5) || (drmu2j1 < 0.5) || (drmu2j2 < 0.5) || (drpho1j1 <  0.5) || (drpho1j2 < 0.5) || (drpho2j1 < 0.5) || (drpho2j2  < 0.5))  continue;
*/
    if(elec1->Charge == elec2 ->Charge) continue;

/*  I wnat my jets to have a DeltaR min of 0.5 between electrons and jets */
    dre1j1 = (elec1->P4()).DeltaR(bjet1->P4());
    dre1j2 = (elec1->P4()).DeltaR(bjet2->P4());
    dre2j1 = (elec2->P4()).DeltaR(bjet1->P4());
    dre2j2 = (elec2->P4()).DeltaR(bjet2->P4());
    drj1j2 = (bjet1->P4()).DeltaR(bjet2->P4());

    if ((dre1j1 < 0.5) || (dre1j2 < 0.5) || (dre2j1 < 0.5) || (dre2j2 < 0.5)) continue;
/*      Ends of the Delta R cuts        */
    MassEE = ((elec1->P4()) +(elec2->P4())).M();


    EnergyEE= ((elec1->P4()) +(elec2->P4())).E();
    MassBB= ((bjet1->P4()) +(bjet2->P4())).M();
    RecoilMass = TMath::Sqrt( EnergyCM*EnergyCM + MassEE*MassEE - 2*EnergyCM*EnergyEE);
    Met = (MissingET*) branchMissingET->At(0);
    //Met = (MissingET *) branchMET->At(0);
    //histMET->Fill(Met->MET);


    if(cuts==0 && ((RecoilMass > 135) || (RecoilMass < 120))) continue; 
    if(cutsee==0 && ((MassEE < 80) || (MassEE > 100))) continue; //s-channel vs t-channel
    else if(cutsee == 1 && EnergyCM == 365 && ( MassEE < 110)) continue;
    //Energy4 = ((elec1->P4()) +(elec2->P4()) + (bjet1->P4()) + (bjet2->P4())+ (muon1->P4()) + (muon2->P4()) + (photon1->P4()) + (photon2->P4())).E();
    Energy4 = ((elec1->P4()) +(elec2->P4()) + (bjet1->P4()) + (bjet2->P4())).E();
    hist1->Fill(RecoilMass,sc);
  }
}

void TesteE4CutsGrace(const Int_t EnergyCM)
{
  const char *inputFile;
  const char *inputFile1;
  const char *inputFile2;
  const char *inputFile3;
  const char *inputFile4;
  const char *inputFile8;

  Int_t Decay;
  Int_t cuts, cutsee;
  cout << "\nType of data:\n\tL.O data[0]\n\tLoop-induced data[1]" << endl;
  Decay=0;

  cout << "Recoil Mass Cuts:\n\tWith cuts:0\n\tWithout cuts:1"<< endl;
  cin >> cuts;

  cout <<"ee Invariant Mass cuts:\n\ts-channel (240 & 365):0\n\tt-channel (365):1\n\tnone:2" << endl;
  cin >> cutsee;

  if (EnergyCM == 240)
  {
    if (Decay == 0)
    {
      const char *x1 = "/root/src/delphes-master/Leading_order240.root";
      const char *x2 = "/root/src/delphes-master/NextLeading_order240.root";
      


      inputFile= x1;
      inputFile1= x2;

    }
          
  }
  else if ( EnergyCM == 365)
  {
    if(Decay == 0)
    {
      const char *y1 = "/root/src/delphes-master/Leading_order365.root";
      const char *y2 = "/root/src/delphes-master/NextLeading_order365.root";

      inputFile= y1;
      inputFile1= y2;

    }
    
  }

  Double_t HistMax,HistMin;
  if(cuts==0) 
  {
    HistMax = 240;
    HistMin = 0;
  }
  else if(cuts==1) 
  {
    HistMax = 200;
    HistMin = 0;
  }
  Double_t RangeX=(HistMax-HistMin);
  Double_t BinperGeV= 3; // in 1 bin I have GeVperBin GeV's
  Double_t GeVperBin=1/BinperGeV;
  Int_t Nbins = RangeX*BinperGeV;

  cout << " Number of Bins: " << Nbins << endl;

  
  TH1F *histE4 = new TH1F("E4","",Nbins,HistMin,HistMax);
  TH1F *histE41 = new TH1F("E41","",Nbins,HistMin,HistMax);


  histE4->Sumw2();
  histE41->Sumw2();
 

/* ------------ Luminosidade -----------------*/
  TFile *Lumi = new TFile("/root/src/delphes-master/Codes/Lumi.root", "READ");
  TTree *tree = (TTree*)Lumi->Get("Luminosidade");

  Double_t Lumi1, Lumi2, LumiXsc;

  tree->SetBranchAddress("eeHee",&Lumi1);
  tree->SetBranchAddress("eeHee_noborn",&Lumi2);
  if (EnergyCM == 240)
  {
    tree->GetEntry(0);

  }
  else if (EnergyCM == 365)
  {
    tree->GetEntry(1);
  }
  cout << endl;
  cout << "Lumi:" << endl;
  cout << "\teeHee LO: "<< Lumi1 << "\n\teeHee NLO:" << Lumi2  << endl;
  
  Lumi->Close();

/* ---------    Scale ---------*/
  Double_t scl, Luminosidade1, scl1, scl2, scl3,scl5,scl6,scl7,scl8;
  if (EnergyCM == 240)
  {
    scl = 5000000.0 / Lumi1;//eeHeetaskc_240_LO
    scl1 = 5000000.0 / Lumi2;//eeHee_240_NLO
    Luminosidade1 = 5;
    LumiXsc=5000000.0;
  }
  else if( EnergyCM == 365 )
  {
    scl = 1500000.0 / Lumi1;//eeHeetaskc_365_LO
    scl1 = 1500000.0 / Lumi2;//eeHee_365_NLO
    Luminosidade1 = 1.5;
    LumiXsc=1500000.0;
  }
/*----------- ------------  ---------*/


gSystem->Load("libDelphes");

//eeHee
TChain chain("Delphes");
chain.Add(inputFile);
ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
AnalyseEvents(treeReader,histE4,EnergyCM,cuts,cutsee,scl);
delete treeReader;

//eeHee_noborn
TChain chain3("Delphes");
chain3.Add(inputFile1); 
ExRootTreeReader *treeReader3 = new ExRootTreeReader(&chain3);
AnalyseEvents(treeReader3,histE41,EnergyCM,cuts,cutsee,scl1);
delete treeReader3;


/*---------------------          CANVAS      ----------------*/
TCanvas *c1 = new TCanvas("c1","Energy 4part", 1200, 1200);

gStyle->SetOptStat(0);
gPad->SetTicks(1,1);
c1->SetGrid();
c1->cd(1);


histE4->SetLineWidth(2);
histE4->SetLineColor(kGreen-1);

histE41->SetLineWidth(2);
histE41->SetLineColor(kPink+3);



histE4->SetTitle("Differential cross section with recoil mass");
histE4->GetXaxis()->SetTitle("M_{recoil} [GeV]");
histE4->GetYaxis()->SetTitle(Form("#frac{d#sigma}{dM_{recoil}} at %.1f ab^{-1}", Luminosidade1));

histE4->Scale(1/(LumiXsc*GeVperBin*Nbins)); // Differential cross-section dXsecton/(FCC_luminosity*GeVperBin*NumberofBins)
histE41->Scale(1/(LumiXsc*GeVperBin*Nbins));

histE4->GetXaxis()->SetRangeUser(115, 165);
histE41->GetXaxis()->SetRangeUser(115, 165);

histE4->Draw("hist");
histE41->Draw("hist same");

Double_t IntegEE = histE4->Integral();
Double_t IntegEE1 = histE41->Integral();

cout <<"\nIntegral eeHee_LO: " << IntegEE;
cout <<"\nIntegral eeHee_NLO: " << IntegEE1;

TLegend *leg1 = new TLegend(0.55,0.75,0.9,0.9);
if(Decay == 0)leg1->SetHeader(Form("#sqrt{s} = %d GeV wo Higgs decay",EnergyCM),"C");


if(cuts==0)leg1->AddEntry((TObject*)0,"Cuts: M_{recoil}#in[120,135] GeV","");
if(cutsee==0)leg1->AddEntry((TObject*)0,"Cuts: M_{ee}#in[80,100] GeV","");
else if(cutsee==1)leg1->AddEntry((TObject*)0,"Cuts: M_{ee} > 110 GeV","");

leg1->AddEntry(histE4, Form("eeHee_LO: Mean=%g, Std Dev=%g", histE4->GetMean(),histE4->GetStdDev()),"l");
leg1->AddEntry(histE41, Form("eeHee_NLO: Mean=%g, Std Dev=%g", histE41->GetMean(),histE41->GetStdDev()),"l");
/*

leg1->AddEntry(histE4, Form("eeHee_LO"),"l");
leg1->AddEntry(histE43, Form("eeHee_NLO"),"l");

*/
leg1->Draw("SAME E");

}

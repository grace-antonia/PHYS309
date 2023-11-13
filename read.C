#include <iostream>
#include <fstream>
#include <string>


void read()
{
    fstream file;
    file.open("Lumi75.txt", ios::in);

    
    Double_t x5, x8;

    TFile *output = new TFile("Lumi75.root", "recreate");

    TTree *tree = new TTree("Luminosidade", "Lumi");

    tree->Branch("eeHee", &x5);
    tree->Branch("eeHee_noborn", &x8);
    

    while(1)
    {
        file >> x5 >> x8;


        

        cout << x5 << " "<< x8 << endl;

        tree->Fill();

        if(file.eof()) break;
    }
    


   
    output->Write();
    output->Close();
    file.close();

}





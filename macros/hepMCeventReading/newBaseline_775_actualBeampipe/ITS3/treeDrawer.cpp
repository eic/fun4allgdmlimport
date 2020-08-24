
#include <iostream>
#include <vector>
#include <sstream>

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <THStack.h>
#include <TGraphErrors.h>
#include <TVector.h>




void treeDrawer() {
	
	TFile * infile = new TFile("out_allSi_B_1.5T_FastTrackingEval.root");
	TTree * analysisTree = (TTree *) infile->Get("tracks");  //Tree 
	analysisTree->Print();
	int nEvents = analysisTree->GetEntries();
	std::cout << "Number of events: " << nEvents << std::endl;
	
	//truth
	int eventNumber = 0;
	int particleID = 0;
	int charge = 0;
	float gpx = 0.0;
	float gpy = 0.0;
	float gpz = 0.0;
	float gvx = 0.0;
	float gvy = 0.0;
	float gvz = 0.0;
	//reconstructed
	float px = 0.0;
	float py = 0.0;
	float pz = 0.0;
	
	analysisTree->SetBranchAddress("event", &eventNumber);
	analysisTree->SetBranchAddress("gflavor", &particleID);
	analysisTree->SetBranchAddress("charge", &charge);
	analysisTree->SetBranchAddress("gpx", &gpx);
	analysisTree->SetBranchAddress("gpy", &gpy);
	analysisTree->SetBranchAddress("gpz", &gpz);
	analysisTree->SetBranchAddress("gvx", &gvx);
	analysisTree->SetBranchAddress("gvy", &gvy);
	analysisTree->SetBranchAddress("gvz", &gvz);
	analysisTree->SetBranchAddress("px", &px);
	analysisTree->SetBranchAddress("py", &py);
	analysisTree->SetBranchAddress("pz", &pz);
	
	int currentEvent = 0; //Keeping track of the current event number
	std::vector<int> particleIDs = {};
	std::vector<double> pionMomenta = {};
	std::vector<std::vector<double>> pionMomComp = {}; //Holding the components
	std::vector<double> kaonMomenta = {};
	std::vector<std::vector<double>> kaonMomComp = {};
	
	std::vector<int> pionCharge = {};
	std::vector<int> kaonCharge = {};
	
	
	double kaonMass = 0.494;
	double pionMass = 0.140;
	
	//Histogram for the D0 invariant mass. Or rather the pion-kaon combinations
	TH1D * invMass = new TH1D("D0 mass", "Pion-kaon invariant mass", 100, 1.7, 2);
	
	
	for (int i = 0; i < nEvents; i++) {
		analysisTree->GetEntry(i);
		//If we have a new event
		if (eventNumber != currentEvent) {
			
			if (kaonMomenta.size() != 0) {
				//for (auto kaonMom : kaonMomenta) {
				for (unsigned kIndex = 0; kIndex < kaonMomenta.size(); kIndex++) {
					double kaonMom = kaonMomenta[kIndex];
					double kaonEnergy = sqrt(kaonMom*kaonMom + kaonMass*kaonMass);
					
					std::vector<double> kComponents = kaonMomComp[kIndex];
					//loop over the pions
					//for (auto pionMom : pionMomenta) {
					for (unsigned pIndex = 0; pIndex < pionMomenta.size(); pIndex++) {
						//If they have the same charge: just skip.
						if (pionCharge[pIndex] == kaonCharge[kIndex]) {
							continue;
						}
						
						double pionMom = pionMomenta[pIndex];
						double pionEnergy = sqrt(pionMom*pionMom + pionMass*pionMass);
						
						std::vector<double> pComponents = pionMomComp[pIndex];
						double momDotProd = pComponents[0]*kComponents[0] + pComponents[1]*kComponents[1] + pComponents[2]*kComponents[2];
						std::cout << "Dot product: " << momDotProd << std::endl;
						double d0MassSquared = pionMass*pionMass + kaonMass*kaonMass + 2*(kaonEnergy*pionEnergy - momDotProd);
						std::cout << "D0 mass: " << sqrt(d0MassSquared) << std::endl;
						if (sqrt(d0MassSquared) > 1.7 && sqrt(d0MassSquared) < 2.0) {
							invMass->Fill(sqrt(d0MassSquared));
						}
					}
					
					//std::cout << mom << " ";
				}
				//std::cout << std::endl;
				//for (auto mom : pionMomenta) {
				//	std::cout << mom << " ";
				//}
				//std::cout << std::endl;
			}
			
			currentEvent = eventNumber;
			std::cout << eventNumber << std::endl;
			
			//Clear up vectors for next event
			particleIDs = {};
			pionMomenta = {};
			kaonMomenta = {};
			pionMomComp = {};
			kaonMomComp = {};
			pionCharge = {};
			kaonCharge = {};
		}
		
		particleIDs.push_back(particleID);
		if (abs(particleID) == 211 && sqrt(gpx*gpx + gpy*gpy) > 0.25) {
			//pionMomenta.push_back(sqrt(gpx*gpx + gpy*gpy + gpz*gpz));
			//pionMomComp.push_back({gpx, gpy, gpz});
			pionMomenta.push_back(sqrt(px*px + py*py + pz*pz));
			pionMomComp.push_back({px, py, pz});
			pionCharge.push_back(charge);
		}
		if (abs(particleID) == 321 && sqrt(gpx*gpx + gpy*gpy) > 0.25) {
			//kaonMomenta.push_back(sqrt(gpx*gpx + gpy*gpy + gpz*gpz));
			//kaonMomComp.push_back({gpx, gpy, gpz});
			kaonMomenta.push_back(sqrt(px*px + py*py + pz*pz));
			kaonMomComp.push_back({px, py, pz});
			kaonCharge.push_back(charge);
		}
		
		
	}
	
	invMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
	invMass->GetYaxis()->SetTitle("Counts");
	invMass->Draw();
	
	//TH2D * baseHistogram = new TH2D("baseHistogram", "baseHistogram", 1000, -0.2, 0.2, 1000, -9, 9);
	////gStyle->SetOptStat(111100);
	//analysisTree->Draw("pcaz:pcax >> baseHistogram", "", "COLZ");
	//baseHistogram->GetXaxis()->SetTitle("pcax");
    //baseHistogram->GetYaxis()->SetTitle("pcaz");
    //
    ////Let's make one with all of them as well, why not
	//TCanvas * c3 = new TCanvas("c3", "c3",900,900);
	//gStyle->SetOptStat(0);
	//
	//TPad * topPad = new TPad("topPad", "topPad", 0.0, 0.5, 1.0, 1.0);
	//TPad * leftPad = new TPad("topPad", "topPad", 0.0, 0.0, 0.5, 0.5);
	//TPad * rightPad = new TPad("topPad", "topPad", 0.5, 0.0, 1.0, 0.5);
	//
	//topPad->Draw();
	//leftPad->Draw();
	//rightPad->Draw();
	//
	//topPad->cd();
	//baseHistogram->Draw("COLZ");
	//
	//leftPad->cd();
	//baseHistogram->ProjectionY()->Draw();
	//
	//rightPad->cd();
	//baseHistogram->ProjectionX()->Draw();
    
}

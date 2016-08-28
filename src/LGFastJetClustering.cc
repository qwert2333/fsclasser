#include "LGFastJetClustering.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "EVENT/LCFloatVec.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <UTIL/LCTOOLS.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

LGFastJetClustering aLGFastJetClustering ;

LGFastJetClustering::LGFastJetClustering() : Processor("LGFastJetClustering") {

	// Processor description
	_description = "FastJet Clustering ..." ;


	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputCollection",
			"Collection of reconstructed particles",
			_inputCollection,
			std::string("Unset") );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputMCTruthMap",
			"Collection of reconstructed particles",
			_inputMCTruthMap,
			std::string("Unset") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"OutputCollection",
			"Name of collection with the found jets",
			_outputCollection,
			std::string("Unset") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"IsoLepCollection",
			"Name of collection with the Isolated Lepton",
			_outputIsoLepCol,
			std::string("IsoLepCol") );
	
	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"RemainCollection",
			"Name of collection without the Isolated Lepton",
			_outputRemainCol,
			std::string("RemainPFOs") );

	registerProcessorParameter("Debug",
			"debug printout",
			_print,
			int(0)); 

	registerProcessorParameter("Algorithm",
			"FastJet algorithm",
			sAlgorithm,
			std::string("antikt_algorithm")); 

	registerProcessorParameter("R",
			"R Parameter",
			_RPar,
			double(0.7)); 
	
	registerProcessorParameter("P",
			"P Parameter",
			_PPar,
			double(0.0)); 

	registerProcessorParameter("EjetMin",
			"Ejet",
			_eJet,
			double(10.0)); 

	registerProcessorParameter("PtMin",
			"PtMin Cut",
			_PtCut,
			double(0.0)); 

	registerProcessorParameter("NJets",
			"max nb of jets",
			_nJetMax,
			int(25)); 

	registerProcessorParameter("nPFOmin",
			"min nb of input PFOs",
			_nPFOmin,
			int(5)); 

	registerProcessorParameter("InclusiveExclusive",
			"Inclusive = 0; Exclusive != 0",
			_InclusiveExclusive,
			int(1)); 

	registerProcessorParameter("FillTree",
			"tuple",
			_fillTree,
			int(0)); 

	registerProcessorParameter("RemoveList",
			"partilces Removed before JetClustering",
			_RemoveList,
			vector<string>()); 
	
	registerProcessorParameter("EnergyCut",
			"Energy Threshold of particles",
			_EnergyCut,
			double(0)); 

}

void LGFastJetClustering::init() { 

	printParameters() ;

	if(_fillTree){

		if(_print>1)cout << "LGFastJetClustering: Making Tuples" << endl;
		_rootfile = new TFile("LGFastJetClustering.root","RECREATE");
		_rootfile->cd("");

		// Declaration of Tree 
		_Etree = new TTree("Events","DST Events");
		_Etree->Branch("NRun",&_nRun,"NRun/I");
		_Etree->Branch("NEvt",&_nEvt,"NEvt/I");
		_Etree->Branch("Ecms",&_eCMS,"Ecms/F");
		_Etree->Branch("NJets",&_nJets,"NJets/I");
		_Etree->Branch("NJetsHE",&_nJetsHE,"NJetsHE/I");
		_Etree->Branch("RPar",&_rp,"RPar/D");
		_Etree->Branch("jetVector",&_jetVector,"jetVector[5][50]/F");

	}

	if(sAlgorithm=="kt_algorithm"){
		SetAlgorithm(fastjet::kt_algorithm);
	}else if(sAlgorithm=="cambridge_algorithm"){ 
		SetAlgorithm(fastjet::cambridge_algorithm);
	}else if(sAlgorithm=="antikt_algorithm"){ 
		SetAlgorithm(fastjet::antikt_algorithm);
	}else if(sAlgorithm=="ee_kt_algorithm"){ 
		SetAlgorithm(fastjet::ee_kt_algorithm);
	}else if(sAlgorithm=="genkt_algorithm"){ 
		SetAlgorithm(fastjet::genkt_algorithm);
	}else if(sAlgorithm=="ee_genkt_algorithm"){ 
		SetAlgorithm(fastjet::ee_genkt_algorithm);
	}else if(sAlgorithm=="plugin_algorithm"){ 
		SetAlgorithm(fastjet::plugin_algorithm);
	}else{
       printf("please chooose the proper algorithm, now it is %s\n", sAlgorithm.c_str());
		 exit(1);
	}

	RemoveResonance = -9999; 

	if ( _RemoveList.size()>0){
		RemoveResonance = atof(_RemoveList[0].c_str());
		_Charged = 0 ;
		for( unsigned int i =1; i< _RemoveList.size();i++)
		{
			_PdgidList.push_back(atoi(_RemoveList[i].c_str()));
			if(_PdgidList[i-1]!=0) _Charged=1;
		}
	}
}

void LGFastJetClustering::processRunHeader( LCRunHeader* run) { 

	_nRun = run->getRunNumber() ;

} 

void LGFastJetClustering::processEvent( LCEvent * evt ) { 

	_jetsCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	_leptCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	_leftCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	_jetsCol->setSubset(true) ;
	_leptCol->setSubset(true) ;
	_leftCol->setSubset(true) ;

	_nRun = evt->getRunNumber();
	_nEvt = evt->getEventNumber();

	if(_print>0)cout <<"Run " << _nRun << " Evt " << _nEvt << endl;

	for(int i1=0;i1<5;i1++){
		for(int i2=0;i2<50;i2++){
			_jetVector[i1][i2]=0.;
		}
	}

	LCCollection* enflowcol  =evt->getCollection(_inputCollection);
	LCCollection* MCTruthMap =evt->getCollection(_inputMCTruthMap);
	LCRelationNavigator *navMCTL   = new LCRelationNavigator(MCTruthMap);
	
	int nenflow =  enflowcol->getNumberOfElements(); 
	if(_print>0)cout <<" # of tracks & clusters " << nenflow << endl;

	double px, py, pz, E;

	fastjet::JetAlgorithm algorithm = GetAlgorithm();
	if(_print>0)cout <<" Jet algorithm is  " << algorithm << endl;
	
	vector<fastjet::PseudoJet> input_particles;

	Double_t dMass=9999; 
	vector<int> plist; plist.clear();
	int pdgid1=0, pdgid2=0;
	int PDGID1=0, PDGID2=0;
	ReconstructedParticle *save_lep1=0, * save_lep2=0; 
	if ( RemoveResonance >0 ) {
		if(_print>0)cout << " RemoveList " << _RemoveList[0] << " "<< _RemoveList[1] << " " << _RemoveList[2] << endl;
		for ( int i=0; i<nenflow-1 ; i++){
			ReconstructedParticle* enflow1 = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( i ));
			Double_t charge1 = enflow1->getCharge();
			if( (int)fabs(charge1) != _Charged) continue;
			//
			LCObjectVec vecMCTL1            = navMCTL->getRelatedToObjects(enflow1);
			if (vecMCTL1.size() > 0) {
				for(unsigned int k=0; k< vecMCTL1.size(); k++){
					MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL1[k]);
					pdgid1=mcp->getPDG();
				}
				if(_print>0)cout << " pdgid1 " << pdgid1 << " "<< vecMCTL1.size() << endl;
			}
			int skip1=1;
			for(unsigned int s =0; s<_PdgidList.size(); s++)
				if( pdgid1 == _PdgidList[s]) skip1=0;
			if(skip1) continue;
			//
			Double_t energy1 = enflow1->getEnergy();
			TVector3 momentum1 = TVector3(enflow1->getMomentum());
			TLorentzVector p41 = TLorentzVector(momentum1,energy1);
			//
			for ( int j=i+1; j<nenflow ; j++){
				ReconstructedParticle* enflow2 = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( j ));
				Double_t charge2 = enflow2->getCharge();
				if( (int)fabs(charge2)!=_Charged) continue;
				if( fabs(charge1+charge2) > 1e-6) continue;
				//
				LCObjectVec vecMCTL2            = navMCTL->getRelatedToObjects(enflow2);
				if (vecMCTL2.size() > 0) {
					for(unsigned int k=0; k< vecMCTL2.size(); k++){
						MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL2[k]);
						pdgid2=mcp->getPDG();
					}
					if(_print>0)cout << " pdgid2 " << pdgid2 << " "<< vecMCTL2.size() << endl;
				}
				int skip2=1;
				for(unsigned int s =0; s<_PdgidList.size(); s++)
					if( pdgid2 == _PdgidList[s]) skip2=0;
				if(skip2) continue;
				//
				Double_t energy2 = enflow2->getEnergy();
				TVector3 momentum2 = TVector3(enflow2->getMomentum());
				TLorentzVector p42 = TLorentzVector(momentum2,energy2);
				double dmass = fabs( (p41+p42).M() - RemoveResonance ); 
				if ( dmass < dMass){
					save_lep1= enflow1;
					save_lep2= enflow2;
					dMass=dmass;
					PDGID1=pdgid1;
					PDGID2=pdgid2;
					plist.clear();
					plist.push_back(i);
					plist.push_back(j);
				}
			}
		}
	}
	if(save_lep1!=0) _leptCol->addElement(save_lep1);
	if(save_lep2!=0) _leptCol->addElement(save_lep2);
	for ( int i=0; i<nenflow ; i++){
		ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( i ));
		if(save_lep1!=enflow && save_lep2!=enflow) _leftCol->addElement(enflow);
	}
	//
	if(_print>1)cout << "PDGID1 & PDGID2 " << PDGID1<<" "<<PDGID2 <<" dMass = "<<dMass<< endl;
	delete navMCTL;
	if(_print>1)cout << "# of RemoveList is " << plist.size() << endl;
	//
	//
	//
	for ( int ienflow=0; ienflow<nenflow ; ienflow++){
		ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( ienflow ));
		if ( enflow->getEnergy()<_EnergyCut ) continue;;
		int skip = 0;
		for( unsigned int i=0; i< plist.size(); i++){
	        if ( ienflow == plist[i] ) skip++;	
		}
		if( skip == 0 ){
			px = enflow->getMomentum()[0];
			py = enflow->getMomentum()[1];
			pz = enflow->getMomentum()[2];
			E  = enflow->getEnergy();

			fastjet::PseudoJet thisPtc(px,py,pz,E);
			thisPtc.set_user_index(ienflow);
			if(_print>2) printf("id, px, py, pz, E= %5d, %10.4f, %10.4f, %10.4f, %10.4f\n", enflow->id(),px, py, pz, E);

			input_particles.push_back(thisPtc);
		}
	}

	if(_print>1) printf("%3d PFOs used for jet-clustering\n", (int)input_particles.size());
	double _ymin[8];
	memset(_ymin, 0, sizeof(double)*8 );
	
	vector<ReconstructedParticleImpl*> trash;
	if( (int)input_particles.size() > _nPFOmin ) {
		_rp = _RPar;
		_pp = _PPar;

		fastjet::Strategy strategy = fastjet::Best;

		fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;

		float momentum[3], energy, mass;

		fastjet::JetDefinition* jet_def =NULL;

		if(algorithm==fastjet::kt_algorithm || algorithm==fastjet::antikt_algorithm || algorithm==fastjet::cambridge_algorithm){
			jet_def = new fastjet::JetDefinition(algorithm, _rp, recomb_scheme, strategy); 
		}else if(algorithm==fastjet::cambridge_algorithm){
			jet_def = new fastjet::JetDefinition(algorithm, _rp, recomb_scheme, strategy);
		}else if(algorithm==fastjet::ee_kt_algorithm){ 
			jet_def = new fastjet::JetDefinition(algorithm); 
		}else if(algorithm==fastjet::genkt_algorithm || algorithm==fastjet::ee_genkt_algorithm ){ 
			jet_def = new fastjet::JetDefinition(algorithm, _rp, _pp, recomb_scheme, strategy); 
		}else if (algorithm==fastjet::plugin_algorithm ){
			double cone_radius = _rp;
			double overlap_threshold = _pp;
			fastjet::SISConePlugin siscone(cone_radius, overlap_threshold);
			jet_def = new fastjet::JetDefinition(& siscone);
		}

		fastjet::ClusterSequence cs(input_particles, *jet_def);

		vector<fastjet::PseudoJet> JetsVec ;

		if( _InclusiveExclusive ) JetsVec = cs.exclusive_jets(_nJetMax);
		else                      JetsVec = cs.inclusive_jets(_PtCut  );
		
		vector<fastjet::PseudoJet> sortedJets = sorted_by_E(JetsVec);

	
		for(int iy=1; iy<9;iy++){	
			_ymin[iy-1] = cs.exclusive_ymerge (iy); // the ymin corresponding to the recombination that went from iy+1 to iy jets
			if(_print>1)printf("%3d,  ymin = %12.5f\n", iy, _ymin[iy-1] );
		}	
		

		int nmx = sortedJets.size();
		if( _InclusiveExclusive>0 && nmx > _nJetMax) nmx = _nJetMax;

		_nJets = sortedJets.size();

		delete jet_def; 


		if(_print>0)printf("LGFastJetClustering: Nb of Jets %3d\n", (int)sortedJets.size() );
		_nJetsHE=0;
		for(int ij=0; ij<_nJets;ij++){

			momentum[0] = sortedJets[ij].px();
			momentum[1] = sortedJets[ij].py();
			momentum[2] = sortedJets[ij].pz();
			energy      = sortedJets[ij].e();
			mass        = sortedJets[ij].m();

			if(_print>1) printf( "Jet %2d  Energy = %10.4f\n", ij,  energy );

			if(energy>_eJet && ij < 10 ){
				
				vector<fastjet::PseudoJet> jetConstituents = cs.constituents(sortedJets[ij]);
				//if( jetConstituents.size() < _nJetMax ) continue;

				_jetVector[0][ij] = sortedJets[ij].px();
				_jetVector[1][ij] = sortedJets[ij].py();
				_jetVector[2][ij] = sortedJets[ij].pz();
				_jetVector[3][ij] = sortedJets[ij].e();
				_jetVector[4][ij] = (float)jetConstituents.size();

				_nJetsHE++;

				if( _InclusiveExclusive==0 || ij<_nJetMax){
					ReconstructedParticleImpl* Jets = new ReconstructedParticleImpl;
					for(unsigned ip = 0; ip < jetConstituents.size(); ip++){
						if(jetConstituents[ip].user_index()>=0 && jetConstituents[ip].user_index() < enflowcol->getNumberOfElements()){
							ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt(jetConstituents[ip].user_index()));
							Jets->addParticle(enflow);
							if(_print>2) cout << "add to particle list  " << enflow->id() << endl;
						}
					}
					Jets->setMomentum(momentum);
					Jets->setEnergy(energy);
					Jets->setMass(mass);
					Jets->setType(4);
					_jetsCol->addElement(Jets);
					trash.push_back(Jets);	
				}
			}
		}

		if(_fillTree){
			_Etree->Fill();
		}    
	}

	_jetsCol->parameters().setValue( "RPar" , (float)_rp );
	_jetsCol->parameters().setValue( "PPar" , (float)_pp );
	_jetsCol->parameters().setValue( "NJets", (float)_nJets );
	for( int iy=1; iy<9; iy++){
		stringstream yname, dname;
		yname << "y"<<iy<<iy+1;
		dname << "d"<<iy<<iy+1;
		_jetsCol->parameters().setValue( yname.str(), (float)_ymin[iy-1] );
	}

	evt->addCollection(_jetsCol ,_outputCollection) ;
	evt->addCollection(_leptCol ,_outputIsoLepCol) ;
	evt->addCollection(_leftCol ,_outputRemainCol) ;
	//FreeDelAll(trash);	
}

void  LGFastJetClustering::check( LCEvent * evt ) { 
}


void LGFastJetClustering::end(){ 
	std::cout << "LGFastJetClustering::end()  " 
		<< " processed " << _nEvt << " events in " << _nRun << " runs "
		<< std::endl ;
	if(_fillTree){
		if(_print>1)cout << "LGFastJetClustering: Saving Tuples" << endl;
		_rootfile->cd("");
		_rootfile->Write();
		_rootfile->Close();
	}
}

#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP 

#include "marlin/FastMCParticleType.h"
#include "marlin/ErrorOfSigma.h"


//--- LCIO headers 


#include <iostream>
#include <cmath>

#include "TVector3.h" 
#include "TLorentzVector.h" 
#include "LGFastMCProcessor.h"
#include "LGParticleFactory.h"
#include "LGTrackSmearer.h"
#include "LGClusterSmearer.h"
#include "cepcplotstyle.h"

using namespace lcio ;


namespace marlin{


	LGFastMCProcessor aLGFastMCProcessor ;


	LGFastMCProcessor::LGFastMCProcessor() : Processor("LGFastMCProcessor"),
	_factory(NULL),
	_nRun(-1),
	_nEvt(-1)
	{

		// modify processor description
		_description = "LGFastMCProcessor creates ReconstrcutedParticles from MCParticles " 
			"according to the resolution given in the steering file." ;


		// register steering parameters: name, description, class-variable, default value

		registerInputCollection( LCIO::MCPARTICLE,
				"InputCollectionName" , 
				"Name of the MCParticle input collection"  ,
				_inputCollectionName ,
				std::string("MCParticle") ) ;


		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"RecoParticleCollectionName" , 
				"Name of the ReconstructedParticles output collection"  ,
				_recoParticleCollectionName ,
				std::string("ReconstructedParticles") ) ;

		registerOutputCollection( LCIO::LCRELATION,
				"MCTruthMappingCollectionName" , 
				"Name of the MCTruthMapping output collection"  ,
				_mcTruthCollectionName ,
				std::string("MCTruthMapping") ) ;


		registerProcessorParameter( "MomentumCut" , 
				"No reconstructed particles are produced for smaller momenta (in [GeV])"  ,
				_momentumCut ,
				float( 0.001 ) ) ;

		FloatVec chResDefault ;
		chResDefault.push_back( 5e-5      ) ;
		chResDefault.push_back( 0.00      ) ;
		chResDefault.push_back( 3.1415927 ) ;

		registerProcessorParameter( "ChargedResolution" , 
				"Resolution of charged particles in polar angle range:  d(1/P)  th_min  th_max"  ,
				_initChargedRes ,
				chResDefault ,
				chResDefault.size() ) ;

		FloatVec gammaResDefault ;
		gammaResDefault.push_back( 0.01      ) ;
		gammaResDefault.push_back( 0.10      ) ;
		gammaResDefault.push_back( 0.00      ) ;
		gammaResDefault.push_back( 3.1415927 ) ;

		registerProcessorParameter( "PhotonResolution" , 
				"Resolution dE/E=A+B/sqrt(E/GeV) of photons in polar angle range: A  B th_min  th_max"  ,
				_initPhotonRes ,
				gammaResDefault ,
				gammaResDefault.size() ) ;

		FloatVec hadronResDefault ;
		hadronResDefault.push_back( 0.04      ) ;
		hadronResDefault.push_back( 0.50      ) ;
		hadronResDefault.push_back( 0.00      ) ;
		hadronResDefault.push_back( 3.1415927 ) ;

		registerProcessorParameter( "NeutralHadronResolution" , 
				"Resolution dE/E=A+B/sqrt(E/GeV) of neutral hadrons in polar angle range: A  B th_min  th_max"  ,
				_initNeutralHadronRes ,
				hadronResDefault ,
				hadronResDefault.size() ) ;

		registerProcessorParameter("MakePlots",     "Make some plots for check",  m_makeplots,       0);
		registerProcessorParameter("Smear",         "modeling the detector res.", m_smear,           1);
		registerProcessorParameter("RejectNeutrino","reject the undetectables  ", m_rejectNeutrino,  0);


	}


	void LGFastMCProcessor::init() { 

		// usually a good idea to
		printParameters() ;

		_nRun = 0 ;
		_nEvt = 0 ;



		_factory = 0 ;
#ifdef MARLIN_CLHEP

		LGParticleFactory* simpleFactory  =  new LGParticleFactory() ; 

		simpleFactory->registerIFourVectorSmearer(  new LGTrackSmearer( _initChargedRes ), CHARGED ) ;
		simpleFactory->registerIFourVectorSmearer(  new LGClusterSmearer( _initPhotonRes ), PHOTON ) ;
		simpleFactory->registerIFourVectorSmearer(  new LGClusterSmearer( _initNeutralHadronRes ), NEUTRAL_HADRON ) ;
		simpleFactory->setMomentumCut( _momentumCut ) ;
		simpleFactory->setSmear( m_smear ) ;
		simpleFactory->setNeutrino( m_rejectNeutrino ) ;

		_factory = simpleFactory ;

		streamlog_out( MESSAGE )  << " LGFastMCProcessor::init() : registering LGParticleFactory " << std::endl ;

#endif // MARLIN_CLHEP

		if( m_makeplots>0 ) { 
			h_Momentum[ 0] = new TH1D("p_gamma"      , "Energy   of #gamma  ", 150,  0.0, 150.0);
			h_Momentum[ 1] = new TH1D("p_electron"   , "Momentum of e^{+}   ", 150,  0.0, 150.0);
			h_Momentum[ 2] = new TH1D("p_positron"   , "Momentum of e^{-}   ", 150,  0.0, 150.0);
			h_Momentum[ 3] = new TH1D("p_muonplus"   , "Momentum of #mu^{+} ", 150,  0.0, 150.0);
			h_Momentum[ 4] = new TH1D("p_muonminus"  , "Momentum of #mu^{-} ", 150,  0.0, 150.0);
			h_Momentum[ 5] = new TH1D("p_pionplus"   , "Momentum of #pi^{+} ", 150,  0.0, 150.0);
			h_Momentum[ 6] = new TH1D("p_pionminus"  , "Momentum of #pi^{-} ", 150,  0.0, 150.0);
			h_Momentum[ 7] = new TH1D("p_kaonplus"   , "Momentum of K^{+}   ", 150,  0.0, 150.0);
			h_Momentum[ 8] = new TH1D("p_kaonminus"  , "Momentum of K^{-}   ", 150,  0.0, 150.0);
			h_Momentum[ 9] = new TH1D("p_proton"     , "Momentum of p       ", 150,  0.0, 150.0);
			h_Momentum[10] = new TH1D("p_antiproton" , "Momentum of #bar{p} ", 150,  0.0, 150.0);
			h_Momentum[11] = new TH1D("p_neutron"    , "Momentum of n       ", 150,  0.0, 150.0);
			h_Momentum[12] = new TH1D("p_antineutron", "Momentum of #bar{n} ", 150,  0.0, 150.0);
			h_Momentum[13] = new TH1D("p_Klong"      , "Momentum of K_{L}   ", 150,  0.0, 150.0);
			h_Momentum[14] = new TH1D("p_pizero"     , "Momentum of #pi^{0} ", 150,  0.0, 150.0);
			//
			h_Mass    [ 0] = new TH1D("M_gammagamma" , "Mass of 2#gamma     ", 350,  0.0, 0.700);
			h_Mass    [ 1] = new TH1D("M_pipi"       , "Mass of 2#pi        ", 350,  0.3, 1.000);
			h_Mass    [ 2] = new TH1D("M_pip"        , "Mass of #pi and p   ", 100,  0.9, 1.200);
			h_Mass    [ 3] = new TH1D("M_KK"         , "Mass of 2Kaon       ", 100,  0.9, 1.100);
			h_Mass    [ 4] = new TH1D("M_Jpsi"       , "Mass of 2lepton     ", 100,  2.9, 3.300);
			h_Mass    [ 5] = new TH1D("M_Upsilon"    , "Mass of 2lepton     ", 100,  9.0, 11.00);
			//
			h_CosTheta[ 0] = new TH1D("ISR_Photon"   , "Cos of Theta        ", 200, -1.0,  1.00);
		}
	}


	void LGFastMCProcessor::processRunHeader( LCRunHeader* run) { 
		_nRun++ ;
		FreeDelAll(_ptrash);
		FreeDelAll(_tracktrash);
		FreeDelAll(_clustertrash);
	} 


	void LGFastMCProcessor::processEvent( LCEvent * evt ) { 
		
		FreeDelAll(_ptrash);
		FreeDelAll(_tracktrash);
		FreeDelAll(_clustertrash);

		const LCCollection* mcpCol = evt->getCollection( _inputCollectionName ) ;

		LCCollectionVec * recVec = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;

		LCRelationNavigator relNav( LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;


		vector<TLorentzVector> PhotonList, PiMinusList, PiPlusList, KPlusList, KMinusList, 
			PList, AntiPList, ElectronList, PositronList, MuPlusList, MuMinusList;

		for(int i=0; i<mcpCol->getNumberOfElements() ; i++){

			MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;

			int status      =  (mcp)->getGeneratorStatus(); 
			int pdgid       =  (mcp)->getPDG();
			int nParents    = ((mcp)->getParents()).size();
			int nDaughters  = ((mcp)->getDaughters()).size();
		   int idGParent    = 0, nGGParents=999;
         if ( nParents>0) {
				idGParent = ((mcp)->getParents()[0])->getPDG();
				nGGParents= ((mcp)->getParents()[0])->getParents().size();
			}	
			TVector3 v(mcp->getMomentum());
			double theta   = v.Theta() ;  
			double pmag    = v.Mag();
			int    pdgcode = mcp->getPDG();

			if( m_makeplots>0 ) {
			   if(pdgcode == 22 &&( nParents==0 ||  nParents>0 && idGParent ==22 &&nGGParents==0)) h_CosTheta[0]->Fill(cos(theta), 1.0);	
				if(pdgcode ==   213 ) h_Momentum[ 5]->Fill(pmag,1.0); 
				if(pdgcode ==  -213 ) h_Momentum[ 6]->Fill(pmag,1.0); 
				if(pdgcode ==   323 ) h_Momentum[ 7]->Fill(pmag,1.0); 
				if(pdgcode ==  -323 ) h_Momentum[ 8]->Fill(pmag,1.0); 
				if(pdgcode ==  2214 ) h_Momentum[ 9]->Fill(pmag,1.0); 
				if(pdgcode == -2214 ) h_Momentum[10]->Fill(pmag,1.0); 
				if(pdgcode ==  2114 ) h_Momentum[11]->Fill(pmag,1.0); 
				if(pdgcode == -2114 ) h_Momentum[12]->Fill(pmag,1.0); 
				if(pdgcode ==   313 ) h_Momentum[13]->Fill(pmag,1.0); 
				if(pdgcode ==   111 ) h_Momentum[14]->Fill(pmag,1.0); 
			}
			if( mcp->getGeneratorStatus() == 1 ) { // stable particles only 
				if( m_makeplots>0 ) { 
					if(pdgcode ==    22 ) h_Momentum[ 0]->Fill(pmag,1.0); 
					if(pdgcode ==    11 ) h_Momentum[ 1]->Fill(pmag,1.0); 
					if(pdgcode ==   -11 ) h_Momentum[ 2]->Fill(pmag,1.0); 
					if(pdgcode ==    13 ) h_Momentum[ 3]->Fill(pmag,1.0); 
					if(pdgcode ==   -13 ) h_Momentum[ 4]->Fill(pmag,1.0); 
				}
				if ( fabs(cos(theta) > 0.9995) ) continue;
				ReconstructedParticle*  rec = 0 ;

				if( _factory != 0 ) 
					rec = _factory->createReconstructedParticle( mcp ) ;

				if( rec != 0 ) {
					recVec->addElement( rec ) ;
					relNav.addRelation( rec , mcp ) ;
					TLorentzVector p4(rec->getMomentum(), rec->getEnergy());
					//
					EVENT::TrackVec::const_iterator it_trk = (rec->getTracks()).begin();
					for(; it_trk != (rec->getTracks()).end() ; it_trk++)
						_tracktrash.push_back( *it_trk);
					EVENT::ClusterVec::const_iterator it_clu = (rec->getClusters()).begin();
					for(; it_clu != (rec->getClusters()).end() ; it_clu++)
						_clustertrash.push_back( *it_clu);
					EVENT::ReconstructedParticleVec::const_iterator it_par = (rec->getParticles()).begin();
					for(; it_par != (rec->getParticles()).end() ; it_par++)
						_ptrash.push_back(*it_par);
					//
					if( m_makeplots>0 ) { 
						if(pdgcode ==   22  ) PhotonList  .push_back(p4) ; 
						if(pdgcode ==   11  ) ElectronList.push_back(p4) ; 
						if(pdgcode ==   13  ) MuMinusList .push_back(p4) ; 
						if(pdgcode ==  -11  ) PositronList.push_back(p4) ; 
						if(pdgcode ==  -13  ) MuPlusList  .push_back(p4) ; 
						if(pdgcode ==   213 ) PiPlusList  .push_back(p4) ; 
						if(pdgcode ==  -213 ) PiMinusList .push_back(p4) ; 
						if(pdgcode ==   323 ) KPlusList   .push_back(p4) ; 
						if(pdgcode ==  -323 ) KMinusList  .push_back(p4) ; 
						if(pdgcode ==  2214 ) PList       .push_back(p4) ; 
						if(pdgcode == -2214 ) AntiPList   .push_back(p4) ; 
					}
				}
			}

		}
		recVec->setDefault   ( true  ) ; // only true, false, false and with track/cluster but without particle
		recVec->setSubset    ( true  ) ; // can make reasonble slcio file with PFO collection/LCRelation !!! 
		recVec->setTransient ( true  ) ;

		evt->addCollection( recVec, _recoParticleCollectionName ) ;
		evt->addCollection( relNav.createLCCollection() , _mcTruthCollectionName ) ;
		//
		if( m_makeplots>0 ) {
         
			//printf("No of Photon is %4d\n", PhotonList.size());
			//printf("No of Kaon+  is %4d\n", KPlusList.size());
			//printf("No of Kaon-  is %4d\n", KMinusList.size());
			if ( PhotonList.size()>0){
				for(unsigned int i=0; i<PhotonList.size()-1; i++){
					for(unsigned int j=i+1; j<PhotonList.size(); j++){
						h_Mass[0]->Fill((PhotonList[i]+PhotonList[j]).M());
					}	 
				}
			}

			for(unsigned int i=0; i<PiPlusList.size(); i++){
				for(unsigned int j=0; j<PiMinusList.size(); j++){
					h_Mass[1]->Fill((PiPlusList[i]+PiMinusList[j]).M());
				}	 
			} 
		
			for(unsigned int i=0; i<PiPlusList.size(); i++){
				for(unsigned int j=0; j<AntiPList.size(); j++){
					h_Mass[2]->Fill((PiPlusList[i]+AntiPList[j]).M());
				}	 
			}

			for(unsigned int i=0; i<PiMinusList.size(); i++){
				for(unsigned int j=0; j<PList.size(); j++){
					h_Mass[2]->Fill((PiMinusList[i]+PList[j]).M());
				}	 
			} 

			for(unsigned int i=0; i<KPlusList.size(); i++){
				for(unsigned int j=0; j<KMinusList.size(); j++){
					h_Mass[3]->Fill((KPlusList[i]+KMinusList[j]).M());
				}	 
			}

			for(unsigned int i=0; i<ElectronList.size(); i++){
				for(unsigned int j=0; j<PositronList.size(); j++){
					h_Mass[4]->Fill((ElectronList[i]+PositronList[j]).M());
					h_Mass[5]->Fill((ElectronList[i]+PositronList[j]).M());
				}	 
			} 

			for(unsigned int i=0; i<MuPlusList.size(); i++){
				for(unsigned int j=0; j<MuMinusList.size(); j++){
					h_Mass[4]->Fill((MuPlusList[i]+MuMinusList[j]).M());
					h_Mass[5]->Fill((MuPlusList[i]+MuMinusList[j]).M());
				}	 
			} 

		}	
		//

		_nEvt ++ ;
	}



	void LGFastMCProcessor::check( LCEvent * evt ) { 
            

	}


	void LGFastMCProcessor::end(){ 

		streamlog_out( MESSAGE4 )  << "LGFastMCProcessor::end()  " << name() 
			<< " processed " << _nEvt << " events in " << _nRun << " runs "
			<< std::endl ;

		if( m_makeplots>0 ) { 
			for( int i=0; i<15; i++){
				h_Momentum[i]->Write();
				//
				TH1D *h1 = new TH1D("h1","",h_Momentum[i]->GetNbinsX(), h_Momentum[i]->GetXaxis()->GetXmin(), h_Momentum[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char*)"P(GeV/c)", (char*)"Entries/1.0GeV/c");
				TH1D *h2 = 0; 
				char filename[256], title[256];
				//
				sprintf(filename,"figs/%s", (h_Momentum[i]->GetName()));
				sprintf(title,"%s", (h_Momentum[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Momentum[i], (char*)h_Momentum[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true, false, false, title
					  );
				//
				sprintf(filename,"figs/%s_log", (h_Momentum[i]->GetName()));
				sprintf(title,"%s", (h_Momentum[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Momentum[i], (char*)h_Momentum[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
					  );
				delete h1;
				delete h_Momentum[i];
			}

			for( int i=0; i<6; i++){
				h_Mass[i]->Write();
				//
				TH1D *h1 = new TH1D("h1","",h_Mass[i]->GetNbinsX(), h_Mass[i]->GetXaxis()->GetXmin(), h_Mass[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char*)"P(GeV/c)", (char*)"Entries/1.0GeV/c");
				TH1D *h2 = 0; 
				char filename[256], title[256];
				//
				sprintf(filename,"figs/%s", (h_Mass[i]->GetName()));
				sprintf(title,"%s", (h_Mass[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Mass[i], (char*)h_Mass[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true, false, false, title
					  );
				//
				sprintf(filename,"figs/%s_log", (h_Mass[i]->GetName()));
				sprintf(title,"%s", (h_Mass[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Mass[i], (char*)h_Mass[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
					  );
				delete h1;
				delete h_Mass[i];
			}
			for( int i=0; i<1; i++){
				h_CosTheta[i]->Write();
				//
				TH1D *h1 = new TH1D("h1","",h_CosTheta[i]->GetNbinsX(), h_CosTheta[i]->GetXaxis()->GetXmin(), h_CosTheta[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char*)"P(GeV/c)", (char*)"Entries/1.0GeV/c");
				TH1D *h2 = 0; 
				char filename[256], title[256];
				//
				sprintf(filename,"figs/%s", (h_CosTheta[i]->GetName()));
				sprintf(title,"%s", (h_CosTheta[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_CosTheta[i], (char*)h_CosTheta[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true, false, false, title
						);
				//
				sprintf(filename,"figs/%s_log", (h_CosTheta[i]->GetName()));
				sprintf(title,"%s", (h_CosTheta[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_CosTheta[i], (char*)h_CosTheta[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
						);
				delete h1;
				delete h_CosTheta[i];
			}

		}
	}
}


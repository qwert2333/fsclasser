#ifndef LGFastJetClustering_h
#define LGFastJetClustering_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include "fastjet/JetDefinition.hh"
#include "Utility.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace Utility;

class LGFastJetClustering : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new LGFastJetClustering ; }
    
  LGFastJetClustering() ;
  
  virtual void init() ;

  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
    
  virtual void check( LCEvent * evt ) ; 

  virtual void end() ;

  fastjet::JetAlgorithm  GetAlgorithm()  const {return fAlgorithm;}

  void SetAlgorithm(fastjet::JetAlgorithm f)  {fAlgorithm = f;}

 protected:

  // root file and tree objects
  TFile * _rootfile;
  TTree * _Etree;

  std::string _inputCollection;
  std::string _inputMCTruthMap;
  std::string _outputCollection;
  std::string _outputIsoLepCol;
  std::string _outputRemainCol;
  
  LCCollectionVec* _jetsCol, *_leptCol, *_leftCol;

  fastjet::JetAlgorithm fAlgorithm; 

  std::string sAlgorithm; 

  // px, py, pz, E, nptc
  float _jetVector[5][50];

  float _eCMS;

  int _nRun, _nEvt, _nJets, _nPFOmin, _nJetsHE;

  int _print, _nJetMax, _fillTree;

  vector<string> _RemoveList;
  vector<int>    _PdgidList;
  int            _Charged, _InclusiveExclusive;
  double RemoveResonance, _EnergyCut; 

  double _RPar, _rp, _eJet, _PPar, _pp, _PtCut;
} ;

#endif

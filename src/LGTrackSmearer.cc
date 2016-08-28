#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP / MARLIN_AIDA

#ifdef MARLIN_CLHEP  // only if CLHEP is available !


#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandEngine.h"


#include <cmath>
#include <cstdlib>

#include "LGTrackSmearer.h"

namespace CLHEP{} 
using namespace CLHEP ;



namespace marlin{


	LGTrackSmearer::LGTrackSmearer(const std::vector<float>& resVec ){

		// copy the resolution vector parameters into a more structured vector 

		_resVec.resize(  resVec.size() / ( sizeof( TrackResolution)  / sizeof(float) )  ) ;  // ==3

		int index = 0 ;

		for( unsigned int i=0 ; i <  _resVec.size() ; i++ ){

			float dPP   =  resVec[ index++ ] ;
			float thMin =  resVec[ index++ ] ;
			float thMax =  resVec[ index++ ] ;

			_resVec[i] = TrackResolution( dPP, thMin, thMax ) ;      
		}
	}


	HepLorentzVector LGTrackSmearer::smearedFourVector( const HepLorentzVector& v, int pdgCode ){


		// find resolution for polar angle
		double theta = v.theta() ;  

		if( theta > M_PI && theta < M_PI*2 )  theta = theta - M_PI; // need to transform to [0,pi] 

		double resolution = -1. ; 

		for( unsigned int i=0 ; i <  _resVec.size() ; i++ ){

			if( theta <= _resVec[i].ThMax  &&  theta > _resVec[i].ThMin ) {
				resolution =  _resVec[i].DPP ;
				break ;
			}
		}
		HepLorentzVector sv( 0., 0. , 0., 0. ) ;

		if( resolution > - 1e-10  ) {

			// do the smearing ....

			double P = v.vect().mag() ;

			double deltaP = RandGauss::shoot( 0.0 , P*P*resolution ) ;

			Hep3Vector n3v( v.vect() )  ;

			n3v.setMag( P + deltaP  ) ;

			//       std::cout << " LGTrackSmearer::smearedFourVector P0: " 
			// 		<< P << " - P1 : " << P + deltaP 
			// 		<< " resolution: " << resolution 
			// 		<< std::endl ;


			// assume perfect electron and muon ID and
			// assign pion mass to everything else

			double mass = PION_MASS ; 

			if( std::abs( pdgCode ) == 11  )  // electron

				mass = ELECTRON_MASS ;

			else if( std::abs( pdgCode ) == 13  )  // muon

				mass = MUON_MASS ;

			sv.setVectM(  n3v  , mass  ) ;

		} 

		return sv ;

	}

}

#endif // MARLIN_CLHEP

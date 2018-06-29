#ifndef PHI_HISTOGRAM_H
#define PHI_HISTOGRAM_H


#include "TreeAnalyzer.h"
#include "XmlRange.h"
#include "RooPlotLib.h"
#include "Reporter.h"
#include "FitConfidence.h"
#include "HistoBins.h"

#include "vendor/loguru.h"

#include "TNamed.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TLatex.h"


#include "FemtoDstFormat/BranchReader.h"
#include "FemtoPairFormat/FemtoPair.h"

#include <algorithm>


class PhiHistogram : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

	float nSigmaDeltaY = 3.0;
	float nSigmaDeltaZ = 3.0;
	float deltaTOF     = 1.0;

	HistoBins ptBins;

public:

	const int DEBUG = 1;
	PhiHistogram() {}
	~PhiHistogram() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		nSigmaDeltaY = config.getFloat( "p.nSigmaDeltaY", 3.0 );
		nSigmaDeltaZ = config.getFloat( "p.nSigmaDeltaZ", 3.0 );
		deltaTOF     = config.getFloat( "p.DeltaTOF", 1.0 );

		ptBins.load( config, "bins.pt" );
	}
protected:

	bool passDeltaY( float pt, float dy ){
		double sigy = -17.6867 + 18.4528*exp(0.637142/pt);
		double marg = 0.;
		if(pt>3.) marg = 0.5;

		if( fabs( dy ) <= ( nSigmaDeltaY + marg ) * sigy ) 
			return true;
		
		return false;
	}

	bool passDeltaZ( float pt, float dz ){
		double sigz = -32.6793 + 32.6034 * exp( 0.444217 / pt );
		double marg = 0.;
		if(pt>3.) 
			marg = 0.5;
		if( fabs( dz ) <= ( nSigmaDeltaZ + marg ) * sigz ) 
			return true;
		return false;
	}

	virtual void analyzeEvent(){
		FemtoPair * pair = this->_fpr.get();

		TLorentzVector lv1, lv2, lv;
		lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.105 );
		lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.105 );
		lv = lv1 + lv2;

		// TRIGGER FLAG
		if ( pair->d1_mTriggerFlag <= 0 )
			return;
		if ( pair->d2_mTriggerFlag <= 0 )
			return;

		if ( false == passDeltaY( pair->d1_mPt, pair->d1_mDeltaY ) || false == passDeltaY( pair->d2_mPt, pair->d2_mDeltaY ) )
			return;
		if ( false == passDeltaZ( pair->d1_mPt, pair->d1_mDeltaZ ) || false == passDeltaZ( pair->d2_mPt, pair->d2_mDeltaZ ) )
			return;

		if ( pair->d1_mNSigmaPion < -2 || pair->d1_mNSigmaPion > 3 )
			return;
		if ( pair->d2_mNSigmaPion < -2 || pair->d2_mNSigmaPion > 3 )
			return;

		// if ( pair->d1_mPt < 1.5 && pair->d2_mPt < 1.5 )
		// 	return;

		if ( pair->d1_mPt < 1.1 || pair->d2_mPt < 1.1 )
			return;

		if ( pair->d1_mDCA > 3.0 || pair->d2_mDCA > 3.0 )
			return;

		if ( fabs(pair->d1_mEta) > 0.5 || fabs(pair->d2_mEta) > 0.5 )
			return;
		if ( fabs(lv.Rapidity()) > 0.5 )
			return;
		if ( fabs( pair->d1_mDeltaTimeOfFlight ) > deltaTOF || fabs( pair->d2_mDeltaTimeOfFlight ) > deltaTOF )
			return;

		////////////////////////////////////////////////////////////////////////
		/// Opposite-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 == pair->mChargeSum ){
			book->fill( "uls", lv.M(), lv.Pt() );

			int ipt = ptBins.findBin( pair->mPt );
			// LOG_F( INFO, "ipt = %d", ipt );
			// book->fill( "uls_mass_" +ts(ipt), lv.M() );
			
		} // 0 == pair->mChargeSum
		else {
			book->fill( "ls", lv.M(), lv.Pt() );
		}


	} //analyzeEvent


	virtual void postMake(){
		TreeAnalyzer::postMake();

		RooPlotLib rpl;
		rpl.link( book );

		

	}




};


#endif


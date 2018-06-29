#ifndef PHI_MVA_HISTOGRAM_H
#define PHI_MVA_HISTOGRAM_H


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


class PhiMvaHistogram : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

	HistoBins ptBins;

	float pid = 1.0;
	string yvar;

public:

	const int DEBUG = 1;
	PhiMvaHistogram() {}
	~PhiMvaHistogram() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		ptBins.load( config, "bins.pt" );

		pid = config.get<float>( "p.pid", 1.0 );
		yvar = config.get<string>( "p.y", "pt" );
	}
protected:

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

		// if ( pair->d1_mPt < 1.5 && pair->d2_mPt < 1.5 )
		// 	return;

		if ( pair->d1_mPt < 1.1 || pair->d2_mPt < 1.1 )
			return;

		// if ( fabs(pair->d1_mEta) > 0.5 || fabs(pair->d2_mEta) > 0.5 )
		// 	return;
		if ( fabs(lv.Rapidity()) > 0.5 )
			return;
		

		float pairPid = sqrt( pow( pair->d1_mPid, 2 ) + pow( pair->d2_mPid, 2 ) );
		if ( pairPid < pid )
			return;

		////////////////////////////////////////////////////////////////////////
		/// Opposite-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 == pair->mChargeSum ){
			if ( "pid" == yvar )
				book->fill( "uls", lv.M(), pairPid );
			else 
				book->fill( "uls", lv.M(), lv.Pt() );
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


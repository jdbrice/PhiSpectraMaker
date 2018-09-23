
#include "calc_mean_pT.C"


void draw_spectra( ){

	TFile *f = new TFile( "systematics/total_ratio_systematics.root" );
	TFile *fTBW = new TFile( "/Users/jdb/bnl/work/dimuonAna/data/Cocktail/MTDMTD-phi_mumu.root" );
	TFile *fy = new TFile( "nim-phi-fit-pp-DNN_N6-DNN136.root" );
	// TFile *fy = new TFile( "nim-phi-womega-fit-pp-DNN_N6-DNN136_1p3.root" );

	TH1 *hnYield = (TH1*)fy->Get( "cyield" );

	TF1 * f1TBW = (TF1*)fTBW->Get( "phi_mumu/phi_pT" );

	f1TBW->SetNpx( 10000 );
	TH1 * hTBW = f1TBW->GetHistogram();

	gStyle->SetOptStat( 0 );

	TH1 * hYieldStat = (TH1*)f->Get( "hYieldStat" );
	TH1 * hYieldSys = (TH1*)f->Get( "hYieldSys" );
	TH1 * hTBWSys = (TH1*)f->Get( "hTBWSys" );



	double x[] = { calc_mean_pT( 2.2, 3.0 ), calc_mean_pT( 3.0, 3.5 ), calc_mean_pT( 3.5, 4.0 ), calc_mean_pT( 4.0, 6.0 ), calc_mean_pT( 6.0, 8.0 ) };
	double y[] = { hnYield->GetBinContent( 1 ), hnYield->GetBinContent( 2 ), hnYield->GetBinContent( 3 ), hnYield->GetBinContent( 4 ), hnYield->GetBinContent( 5 ) };
	double exl[] = { calc_mean_pT( 2.2, 3.0 ) - 2.2, calc_mean_pT( 3.0, 3.5 ) - 3.0, calc_mean_pT( 3.5, 4.0 ) - 3.5, calc_mean_pT( 4.0, 6.0 ) - 4.0, calc_mean_pT( 6.0, 8.0 ) - 6.0 };
	double exh[] = { 3.0 - calc_mean_pT( 2.2, 3.0 ), 3.5 - calc_mean_pT( 3.0, 3.5 ), 4.0 - calc_mean_pT( 3.5, 4.0 ), 6.0 - calc_mean_pT( 4.0, 6.0 ), 8.0 - calc_mean_pT( 6.0, 8.0 ) };
	double eyl[] = { hnYield->GetBinError( 1 ), hnYield->GetBinError( 2 ), hnYield->GetBinError( 3 ), hnYield->GetBinError( 4 ), hnYield->GetBinError( 5 ) };
	double eyh[] = { hnYield->GetBinError( 1 ), hnYield->GetBinError( 2 ), hnYield->GetBinError( 3 ), hnYield->GetBinError( 4 ), hnYield->GetBinError( 5 ) };

	double eysl[] = { hYieldSys->GetBinError( 1 ), hYieldSys->GetBinError( 2 ), hYieldSys->GetBinError( 3 ), hYieldSys->GetBinError( 4 ), (hYieldSys->GetBinError( 4 ) / hYieldSys->GetBinContent(4)) * hnYield->GetBinContent(5)  };
	double eysh[] = { hYieldSys->GetBinError( 1 ), hYieldSys->GetBinError( 2 ), hYieldSys->GetBinError( 3 ), hYieldSys->GetBinError( 4 ), (hYieldSys->GetBinError( 4 ) / hYieldSys->GetBinContent(4)) * hnYield->GetBinContent(5)  };


	double  mtd_resp_sys[] = { (0.9-0.6021) * y[0], (0.9-.7783) * y[1], (0.9 - 0.8911) * y[2], 0, 0.05 * y[4] };

	for ( int i = 0; i < 5; i++ ){
		eysl[i] = sqrt( eysl[i]*eysl[i] + mtd_resp_sys[i]*mtd_resp_sys[i] );
		eysh[i] = sqrt( eysh[i]*eysh[i] + mtd_resp_sys[i]*mtd_resp_sys[i] );
	}

	TGraphAsymmErrors *tgestat = new TGraphAsymmErrors( 5, x, y, exl, exh, eyl, eyh );
	TGraphAsymmErrors *tgesys = new TGraphAsymmErrors( 5, x, y, exl, exh, eysl, eysh );


	double Itarget = hTBWSys->Integral( hTBWSys->GetXaxis()->FindBin(2.2), hTBWSys->GetXaxis()->FindBin(3.0) ) * (hTBWSys->GetBinWidth(1) / hTBW->GetBinWidth(1)) ;
	hTBW->Scale( Itarget  / (hTBW->Integral( hTBW->GetXaxis()->FindBin(2.2), hTBW->GetXaxis()->FindBin(3.0) ) ) );
	for ( int i = 1; i <= hTBW->GetXaxis()->GetNbins(); i++ ){
		hTBW->SetBinError( i, hTBW->GetBinContent(i) * 0.26 );
	}

	TCanvas *can = new TCanvas( "can", "", 1600, 900 );
	can->SetTopMargin( 0.03 ); 
	can->SetRightMargin( 0.03 );
	can->SetBottomMargin( 0.13 );
	can->SetLeftMargin( 0.13 );

	hTBW->GetYaxis()->SetTitle( "BR #times #frac{d^{2}N}{dp_{T}dyp_{T} } #frac{1}{2#pi N_{events} }" );
	hTBW->GetXaxis()->SetTitle( "p_{T} (GeV/c)" );
	hTBW->SetMaximum( 1e-6 );
	hTBW->SetMinimum( 1e-12 );
	hTBW->GetXaxis()->SetRangeUser( 1.5, 8.0 );
	hTBW->SetLineColor(kBlue);
	hTBW->SetFillColorAlpha( kBlue, 0.4 );
	hTBW->Draw("pE2");

	tgestat->SetMarkerStyle( 29 );
	tgestat->SetMarkerSize( 3 );
	tgestat->SetMarkerColor( kRed );
	tgestat->Draw( "same p" );

	tgesys->SetMarkerStyle( 29 );
	tgesys->SetMarkerSize( 0.1 );
	tgesys->SetFillColorAlpha( kRed, 0.6 );
	tgesys->Draw( "same E2" );
	can->SetLogy(1);

	TLegend *leg= new TLegend( 0.5, 0.6, 0.97, 0.97 );
	leg->AddEntry( tgestat, "Statistical Uncertainty", "pe" );
	leg->AddEntry( tgesys, "Systematic Uncertainty", "fpe" );
	leg->AddEntry( hTBW, "TBW", "l" );
	leg->AddEntry( hTBW, "TBW dN/dy Uncertainty", "E4" );

	leg->Draw("same");

	can->Print( "present_spectra.pdf" );

	can->SetLogy(1);
	double bins[] = {2.2, 3.0, 3.5, 4.0, 6.0, 8.0};
	TH1 * hratioFrame = new TH1F( "hratioFrame", ";p_{T} (GeV/c); Yield Ratio", 500, 1.5, 8.0 );
	
	for ( int i = 20; i < 30; i++ ){
		hratioFrame->SetBinContent( i, 1.0);
		hratioFrame->SetBinError( i, 0.26 );
	}
	
	

	double ry[5];
	double rey[5];
	double reys[5];


	for ( int i = 0; i < 5; i++ ){
		ry[i] = y[i] / hTBW->GetBinContent( hTBW->GetXaxis()->FindBin( x[i] ) ) ;
		rey[i] = eyh[i] / hTBW->GetBinContent( hTBW->GetXaxis()->FindBin( x[i] ) ) ;
		reys[i] = eysl[i] / hTBW->GetBinContent( hTBW->GetXaxis()->FindBin( x[i] ) ) ;
	}

	TGraphAsymmErrors *tgestatRatio = new TGraphAsymmErrors( 5, x, ry, exl, exh, rey, rey );
	TGraphAsymmErrors *tgesysRatio = new TGraphAsymmErrors( 5, x, ry, exl, exh, reys, reys );


	hratioFrame->SetMinimum( 0.05 );
	hratioFrame->SetMaximum( 10 );
	hratioFrame->SetFillColorAlpha(kGray+2, 0.9);
	hratioFrame->SetLineColor(kGray+2);
	hratioFrame->Draw( "E2");

	tgestatRatio->Draw("same ")	;

	TLine *l = new TLine( 1.5, 1.0, 8.0, 1.0 );
	l->SetLineColor( kBlack );
	l->SetLineStyle( 2 );
	l->SetLineWidth( 2 );
	l->Draw("Same");
	
	tgesysRatio->SetMarkerStyle( 29 );
	tgesysRatio->SetMarkerSize( 0.1 );
	tgesysRatio->SetFillColorAlpha( kRed, 0.6 );
	tgesysRatio->Draw( "same E2" );

	TLegend *leg2= new TLegend( 0.5, 0.8, 0.97, 0.97 );
	leg2->AddEntry( tgestatRatio, "Statistical Uncertainty", "pe" );
	leg2->AddEntry( tgesys, "Systematic Uncertainty", "fpe" );
	leg2->AddEntry( hratioFrame, "TBW dN/dy Uncertainty", "f" );

	leg2->Draw("same");

	can->Print( "present_spectra_ratio.pdf" );


	

}
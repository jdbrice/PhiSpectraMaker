#ifndef PHI_FITTER_H
#define PHI_FITTER_H

#include "HistoAnalyzer.h"
#include "FitConfidence.h"
#include "XmlHistogram.h"

#include "TRandom3.h"
#include "TLatex.h"
#include "TF1.h"

TF1 * phi_fpol = nullptr;
TF1 * phi_fg1 = nullptr;
TF1 * phi_fg2 = nullptr;

Double_t phi_evaluate(Double_t *x, Double_t *p){

	phi_fg1->SetParameters( p[0], p[1], p[2] );
	phi_fpol->SetParameters( p[3], p[4], p[5], p[6], p[7], p[8] );

	return phi_fpol->Eval( x[0] ) + phi_fg1->Eval( x[0] );
}

Double_t phi_and_omega_evaluate(Double_t *x, Double_t *p){

	phi_fg1->SetParameters( p[0], p[1], p[2] );
	phi_fg2->SetParameters( p[3], p[4], p[5] );
	phi_fpol->SetParameters( p[6], p[7], p[8], p[9], p[10], p[11] );

	return phi_fpol->Eval( x[0] ) + phi_fg1->Eval( x[0] ) + phi_fg2->Eval( x[0] );
}


class PhiFitter : public HistoAnalyzer {
protected:
	TRandom3 r;
	int npol;

	XmlHistogram xhMtdR, xhGen, xhTBW;
	TF1 * fitfunction = nullptr;
	
public:


	

	virtual void initialize(){
		HistoAnalyzer::initialize();

		npol = config.get<int>( "p.bgNPol", 4 );
		LOG_F( INFO, "Background: pol%d", npol );
		phi_fpol = new TF1( "phi_fpol", TString::Format( "pol%d", npol ) );
		phi_fg1 = new TF1( "phi_fg1", "gaus" );
		phi_fg2 = new TF1( "phi_fg2", "gaus" );

		xhMtdR.load( config, config.q( nodePath + ".XmlHistogram{name==mtdr_pt}" ) );
		xhGen.load( config, config.q( nodePath + ".XmlHistogram{name==gen_pt}" ) );
		xhTBW.load( config, config.q( nodePath + ".XmlHistogram{name==TBW}" ) );
		

		r.SetSeed(0);
	}



	virtual void phiFit(){
		RooPlotLib rpl;
		rpl.link( book );
		rpl.link( &config );
		gStyle->SetOptStat(0);

		string resoStr = config.get<string>( "p.reso", "phi" );
		string reso2Str = config.get<string>( "p.reso2", "" );
		LOG_F( INFO, "reso2Str=%s", reso2Str.c_str() );
		const char * reso = resoStr.c_str();

		const char * name = config.get<string>( "name", "" ).c_str();
		const char * title = config.get<string>( "title", "" ).c_str();

		book->cd();

		TH2 * huls = get<TH2>( "uls" );
		TH2 * hls = get<TH2>( "ls" );
		LOG_F( INFO, "mass vs pT histo : %p", huls );
		float pt1 = config.get<float>( "p.pt1", 0 );
		float pt2 = config.get<float>( "p.pt2", 1000 );

		Reporter rp( TString::Format( config.get<string>( nodePath + ".output.FitReport:url" ).c_str(), reso, name, pt1, pt2).Data(), 800, 800 );
		rp.margins( 0.05, 0.05, 0.13, 0.15 );
		rp.newPage();
		

		int ipt1 = huls->GetYaxis()->FindBin( pt1 );
		int ipt2 = huls->GetYaxis()->FindBin( pt2 );

		LOG_F( INFO, "%s pT ( %0.3f -> %0.3f ), bins ( %d -> %d )", reso, pt1, pt2, ipt1, ipt2 );

		TH1 * hmass   = huls->ProjectionX( TString::Format( "mass_%d_to_%d", ipt1, ipt2 ), ipt1, ipt2 );
		TH1 * hmassls = hls ->ProjectionX( TString::Format( "massls_%d_to_%d", ipt1, ipt2 ), ipt1, ipt2 );

		hmass->Sumw2();
		hmassls->Sumw2();


		HistoBins sigMassBins;
		sigMassBins.load( config, "bins.mass" );
		TH1 * hmassrb = hmass->Rebin( 
										sigMassBins.nBins(), 
										TString::Format( "mass_%d_to_%d_rb", ipt1, ipt2 ),
										sigMassBins.bins.data() 
									);

		// hmassrb->Scale( 1.0, "width" );
		double bw = hmassrb->GetBinWidth( 1 );

		TH1 * hmassrbls = hmassls->Rebin( 
										sigMassBins.nBins(), 
										TString::Format( "massls_%d_to_%d_rb", ipt1, ipt2 ),
										sigMassBins.bins.data() 
									);

		// hmassrbls->Scale( 1.0, "width" );


		if ( config.get<bool>( "p.rmls", false ) == true ){
			hmassrb->Add( hmassrbls, -1 );
		}

		float mu = 1.013;
		if ( "omega" == resoStr ){
			mu = 0.78;
		} else if ( "k0s" == resoStr ){
			mu = 0.47;
		}

		TF1 * f = nullptr;
		if ( nullptr == fitfunction && "" == reso2Str ){
			LOG_F( INFO, "PHI only" );
			fitfunction = new TF1( "phif", phi_evaluate, 2, 5, 3 + npol );
			fitfunction->SetParameters( 610, mu, 0.015, 1, -1, -1, 3, 2, -1 );

			fitfunction->SetNpx( 500 );

			fitfunction->SetParLimits( 0, 0, 1e9 );
			fitfunction->SetParLimits( 1, mu - 0.1, mu + 0.1);
			fitfunction->SetParLimits( 2, 0.001, 0.06 );
		} else if ( nullptr == fitfunction && "" != reso2Str ){
			LOG_F( INFO, "PHI + OMEGA" );
			fitfunction = new TF1( "phif", phi_and_omega_evaluate, 2, 5, 3 + 3 + npol );
			fitfunction->SetParameters( 310, mu, 0.015, 610, 0.777, 0.017, 1, -1, -1, 3, 2 );
			fitfunction->SetParameter( 11, -1 );

			fitfunction->SetNpx( 500 );

			fitfunction->SetParLimits( 0, 1, 1e3 );
			fitfunction->SetParLimits( 1, mu - 0.051, mu + 0.051);
			fitfunction->SetParLimits( 2, 0.001, 0.06 );

			// fitfunction->SetParLimits( 3, 0, 1e9 );
			fitfunction->SetParLimits( 4, 0.78 - 0.03, 0.78 + 0.01);
			fitfunction->SetParLimits( 5, 0.001, 0.02 );

			fitfunction->FixParameter( 4, 0.777 );
			fitfunction->FixParameter( 5, 0.017 );

			if ( "phi" == resoStr && config.get<bool>( "p.fixPhi" ) ){
				LOG_F( INFO, "FIXING PHI params" );
				fitfunction->FixParameter( 1, 1.016 );
				fitfunction->FixParameter( 2, 0.014 );
			}
		}

		f = fitfunction;
		 // new TF1( "phif", phi_evaluate, 2, 5, 3 + npol );
		rpl.style( f ).set( config, "style.fit" );
		f->SetLineColor(kBlack);
		
		

		// f->SetParameters( 1.2e5, -1.2e5, -1.2e4, 3.7e4, 2.5e4, -1.9e4, 610, mu, 0.015 );
		// f->SetParameters( 1, -1, -1, 3, 2, -1, 610, mu, 0.015 );
		

		

		
		
		float fmin = config.get<float>("p.fitRange:min", 0.85);
		float fmax = config.get<float>("p.fitRange:max", 1.5);;

		if ( "omega" == resoStr ){
			fmin = 0.55;
			fmax = 0.9;
		} else if ( "k0s" == resoStr ){
			fmin = 0.2;
			fmax = 0.9;
		}
		// if ( pt1 < 1.0 ){
		// 	fmax = 4.5;
		// 	fmin = 2.6;
		// }


		
		// hmassrb->Fit( f, "QR", "", fmin, fmax );
		// hmassrb->Fit( f, "QRL", "", fmin, fmax );
		hmassrb->Fit( f, "QRL", "", fmin, fmax );
		TFitResultPtr fresult = hmassrb->Fit( f, "SRLE", "", fmin, fmax );
		
		float fmu  = f->GetParameter( 1 );
		float fsig = f->GetParameter( 2 );
		double Nse = f->IntegralError( fmu - 5*fsig, fmu + 5*fsig ) / bw ;
		if ( 0.0f == Nse ){
			hmassrb->Fit( f, "QR", "", fmin, fmax );
			hmassrb->Fit( f, "R", "", fmin, fmax );
			fresult = hmassrb->Fit( f, "SRLE", "", fmin, fmax );
		}

		fmu  = f->GetParameter( 1 );
		fsig = f->GetParameter( 2 );
		Nse = f->IntegralError( fmu - 3*fsig, fmu + 3*fsig ) / bw;


		float chi2 = f->GetChisquare();
		float ndf = (float)f->GetNDF();
		float rchi2 = chi2 / ndf;
		LOG_F( INFO, "%f / %f = %f", chi2, ndf, rchi2 );
		

		rpl.style( hmassrb ).set( config, "style.mass" ).draw();


		f->SetRange( fmin, fmax );
		phi_fpol->SetRange( fmin, fmax );
		phi_fg1->SetRange( fmin, fmax );

		// TH1 * g = jdb::FitConfidence::fitCL( phi_fg1, "fit_uncert", 0.95, 500 );
		// g->SetLineColor( kRed );
		// g->Draw( "same" );

		rpl.style( f ).set( config, "style.fit" );
		f->SetNpx( 500 );
		
		f->Draw("same");

		rpl.style( phi_fg1 ).set( config, "style.fit" );
		phi_fg1->SetNpx( 500 );
		
		phi_fg1->Draw("same");


		if ( "" != reso2Str ){
			rpl.style( phi_fg2 ).set( config, "style.fit2" );
			phi_fg2->SetNpx( 500 );
			
			phi_fg2->Draw("same");
		}

		phi_fpol->SetLineColor(kRed);
		rpl.style( phi_fpol ).set( config, "style.fitbg" );
		
		phi_fpol->SetNpx( 500 );
		phi_fpol->DrawClone("same");

		// rpl.style( phi_fpol ).set( "style.fitext" );
		// phi_fpol->SetRange( 0.2, 2.5 );
		// phi_fpol->SetNpx( 500 );
		// phi_fpol->Draw("same");

		// TGraphErrors *g = jdb::FitConfidence::choleskyBands( fresult, f );

		

		// float bw = hmassrb->GetBinWidth( 5 );
		cout << "Ns FIT: " << f->GetParameter( 0 ) << endl;
		double Ns = f->Integral( fmu - 3*fsig, fmu + 3*fsig ) / bw;
		double Nbg = phi_fpol->Integral( fmu - 3*fsig, fmu + 3*fsig ) / bw;

		
		double Nbge = phi_fpol->IntegralError( fmu - 3*fsig, fmu + 3*fsig ) / bw;

		// Ns *= 1.0/ hmassrb->GetBinWidth(4);
		// Nbg *= 1.0/ hmassrb->GetBinWidth(4);

		LOG_F( INFO, "Ns=%f", Ns );
		LOG_F( INFO, "Nbg=%f", Nbg );


		Ns = Ns - Nbg;
		Nse = Nse;

		double SoBe = (Ns/Nbg) * sqrt( pow( sqrt(Ns) / Ns,2) + pow( sqrt(Nbg) / Nbg, 2 ) );

		LOG_F( INFO, "S/B=%0.3e #pm %0.3e", Ns/Nbg, SoBe );

		LOG_F( INFO, "Ns-Nbg=%f", Ns );
		
		double N = Ns + Nbg;
		double Ne = sqrt( N );
		double sig = Ns / Ne;
		double sige = sig * (Nse / Ns);

		TLatex lx;
		lx.SetTextSize( 12.0 / 380.0);
		
		lx.DrawLatexNDC( .42, 0.90, "Run15 p+p at #sqrt{s}=200 GeV" );
		
		lx.DrawLatexNDC( .42, 0.85, TString::Format("N_{#phi}^{raw}=%0.3f #pm %0.3f", Ns, Nse) );
		lx.DrawLatexNDC( .42, 0.80, TString::Format("mass=%0.3f #pm %0.3f GeV/c^{2}", f->GetParameter(1), f->GetParError(1)) );
		lx.DrawLatexNDC( .42, 0.75, TString::Format("width=%0.3f #pm %0.3f GeV/c^{2}", f->GetParameter(2), f->GetParError(2)) );
		lx.DrawLatexNDC( .42, 0.70, TString::Format("S/B=%0.3f", Ns/Nbg) );
		lx.DrawLatexNDC( .42, 0.65, TString::Format("S/#sqrt{S + B}=%0.3f", sig) );

		lx.DrawLatexNDC( .72, 0.70, TString::Format("#chi^{2}/NDF=") );

		lx.DrawLatexNDC( .72, 0.65, TString::Format("%0.1f/%0.1f=%0.2f", chi2, ndf, rchi2) );

		// lx.DrawLatexNDC( .18, 0.40, TString::Format("p0=%0.3f #pm %0.2e", f->GetParameter(6), f->GetParError(6) ) );
		
		lx.DrawLatexNDC( .23, 0.28, "|y^{#mu#mu}|<0.5" );
		lx.DrawLatexNDC( .23, 0.23, "|#eta^{#mu}|<0.5 && p^{#mu}_{T} > 1.1 GeV/c" );
		lx.DrawLatexNDC( .23, 0.18, TString::Format( "%s", title ) );

		if ( config.get<bool>("p.showSlice") ){
			lx.DrawLatexNDC( .23, 0.18, TString::Format("%0.2f < %s < %0.2f %s", pt1, config.get<string>("p.y-title", "p_{T}^{#mu#mu}").c_str(), pt2, config.get<string>("p.y-units", "GeV/c").c_str() ) );
		}

		TLegend *leg = new TLegend( 0.23, 0.35, 0.5, 0.55 );
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry( hmassrb, "Data", "lpe" );
		leg->AddEntry( phi_fg1, "Signal" );
		leg->AddEntry( phi_fpol, "Background" );
		if ( config.get<bool>( "p.showLegend", true ) )
			leg->Draw("same");

		rp.savePage();

		float mpt = (pt1 + pt2)/2.0;
		if ( true == config.get<bool>( "p.stepMinOnly", false ) )
			mpt = pt1 + 1e-5;
		if ( book->get( "yield" ) ){
			int ibin = book->get("yield")->GetXaxis()->FindBin( mpt );
			book->setBin( "yield", ibin, Ns, sqrt(Ns) );
		}
		if ( book->get( "mass" ) ){
			int ibin = book->get("mass")->GetXaxis()->FindBin( mpt );
			book->get("mass")->SetBinContent( ibin, f->GetParameter( 1 ) );
			book->get("mass")->SetBinError( ibin, f->GetParError( 1 ) );
			
		}
		if ( book->get( "width" ) ){
			int ibin = book->get("width")->GetXaxis()->FindBin( mpt );
			book->get("width")->SetBinContent( ibin, f->GetParameter( 2 ) );
			book->get("width")->SetBinError( ibin, f->GetParError( 2 ) );
			
		}
		if ( book->get( "SoverB" ) ){
			int ibin = book->get("SoverB")->GetXaxis()->FindBin( mpt );
			book->get("SoverB")->SetBinContent( ibin, Ns/Nbg );
			book->get("SoverB")->SetBinError( ibin, SoBe );
			
		}
		if ( book->get( "significance" ) ){
			int ibin = book->get("significance")->GetXaxis()->FindBin( mpt );
			book->get("significance")->SetBinContent( ibin, sig );
			book->get("significance")->SetBinError( ibin, sige );
			
		}


		f->Write();
	}

	virtual void make(){
		LOG_F( INFO, "" );

		book->cd();
		book->makeAll( nodePath + ".histograms" );

		for ( string n : {"mass", "width", "yield", "SoverB", "significance"} ){
			if ( book->get(n) )
				book->get(n)->Sumw2();
		}

		HistoBins fit_bins( config, "bins.fit" );
		LOG_F( INFO, "Found %d slices to fit", fit_bins.nBins() );
		for ( size_t i = 0; i < fit_bins.nBins(); i++ ){
			LOG_F( INFO, "fit %lu", i );
			float pt1 = fit_bins.getBins()[i];
			float pt2 = fit_bins.getBins()[i+1];
			LOG_F( INFO, "Setting fit y=(%f, %f)", pt1, pt2 );

			config.set( "p.pt1", ts(pt1) );
			if ( false == config.get<bool>( "p.stepMinOnly", false ) )
				config.set( "p.pt2", ts(pt2) );

			phiFit();
		}
	}

	virtual void postMake(){

		TNamed n ( "config", "" );
		n.SetTitle( config.toXml().c_str() );
		n.Write();
		
		Reporter rp( config.get<string>(nodePath + ".output.YieldSummary:url"), 1600, 800 );
		rp.margins( 0.03, 0.05, 0.15, 0.15 );
		rp.newPage();

		RooPlotLib rpl;
		rpl.link( book );
		rpl.link( &config );

		if ( book->get( "yield" ) ){
			// raw yield from fits
			LOG_F( INFO, "e(4) = %f vs. %f", book->get("yield")->GetBinError(4) / book->get("yield")->GetBinContent(4), sqrt(book->get("yield")->GetBinContent(4)) / book->get("yield")->GetBinContent(4)  );
			rpl.style( "yield" ).set( "style.yield" ).draw();
			gPad->SetLogy(1);
			rp.next();

			if ( true == config.get<bool>( "p.scalebw", true ) )
				book->get( "yield" )->Scale( 1.0 / ( 2 * 3.1415926 ), "width" );
			LOG_F( INFO, "e(4) = %f", book->get("yield")->GetBinError(4) / book->get("yield")->GetBinContent(4) );
			if ( true == config.get<bool>( "p.scalept", true ) ){
				for ( int i = 1; i <= book->get("yield")->GetXaxis()->GetNbins(); i++ ){
					float v = book->get("yield")->GetBinContent(i);
					float e = book->get("yield")->GetBinError(i);
					float mpt = book->get("yield")->GetBinCenter(i);
					book->get("yield")->SetBinContent( i, v / mpt );
					book->get("yield")->SetBinError( i, e / mpt );
				}
			}
			LOG_F( INFO, "e(4) = %f", book->get("yield")->GetBinError(4) / book->get("yield")->GetBinContent(4) );
			rpl.style( "yield" ).set( "style.yield" ).draw();
			gPad->SetLogy(1);
			rp.next();
		}


		// load efficiency tables
		// compute total efficiency x acceptance and apply correction
		
		HistoBins ptBins( config, "bins.fit" );
		
		TH1 * hEffAcc = (TH1*)xhMtdR.getTH1().get()->Clone( "hEffAcc" );
		hEffAcc->Divide( xhGen.getTH1().get() );
		hEffAcc->Write();
		TH1 * hrbmtdr_pt = HistoBins::rebin1D( "rbmtdr_pt", xhMtdR.getTH1().get(), ptBins );
		TH1 * hrbgen_pt = HistoBins::rebin1D( "rbgen_pt", xhGen.getTH1().get(), ptBins );

		TH1 * hrbEffAcc = (TH1*)hrbmtdr_pt->Clone( "hrbEffAcc" );
		hrbEffAcc->Divide( hrbgen_pt );
		hrbEffAcc->Write();

		hEffAcc->SetMinimum( 1e-4 );
		hEffAcc->Draw();
		hrbEffAcc->Draw("same");
		rp.next();

		book->clone( "yield", "cyield" );
		TH1 * cyield = book->get("cyield");
		for ( int i = 1; i <= cyield->GetXaxis()->GetNbins(); i++ ){
			float v = cyield->GetBinContent( i );
			float e = cyield->GetBinError( i );
			float eff = hrbEffAcc->GetBinContent( hrbEffAcc->GetXaxis()->FindBin( cyield->GetBinCenter( i ) ) );

			float relEffShift = config.get<float>( "p.RelativeEffShift", 1.0 );
			eff *= relEffShift;

			LOG_F( INFO, "eff = %f at %f", eff, cyield->GetBinCenter( i ) );
			cyield->SetBinContent( i, v / eff );
			cyield->SetBinError( i, (e / eff) );

		}

		
		TH1 * hTBW = xhTBW.getTH1().get();
		hTBW->Scale( (1.73e-2 / hTBW->Integral()) * 2.87e-4, "width" );
		
		
		// hTBW->SetMinimum( 1e-1 );
		rpl.style( hTBW ).set( config, "style.TBW" ).draw();
		hTBW->Write();
		LOG_F( INFO, "e(4) = %f", cyield->GetBinError(4) / cyield->GetBinContent(4) );

		cyield->Scale( 1.0/ ( 3.56e11 ) );
		cyield->Scale( (1.0/0.97) * ( 1.0 / 0.84 ) * (1.0 / 0.55) ); // missing efficiencys, will be formally added to eff
		// i think these were trigger, trigger unit, and the PID efficiency

		LOG_F( INFO, "e(4) = %f", cyield->GetBinError(4) / cyield->GetBinContent(4) );
		

		TH1 * cyieldSys = (TH1*)cyield->Clone( "cyieldSys" );
		for ( int i = 1; i <= cyield->GetXaxis()->GetNbins(); i++ ){
			cyieldSys->SetBinError( i, cyieldSys->GetBinError( i ) * 5 );
		}
		rpl.style( cyieldSys ).set( config, "style.cyieldSys" ).draw();
		rpl.style( cyield ).set( config, "style.cyield" ).draw();

		gPad->SetLogy(1);
		rp.next();

		if ( config.get<float>( "p.RelativeEffShift", 1.0 ) != 1.0 ){
			LOG_F( WARNING, "Shifting Eff by eff*=%f", config.get<float>( "p.RelativeEffShift", 1.0 ) );
		}

		book->clone( "cyield", "yratio" );
		TH1 * yratio = book->get("yratio");
		for ( int i = 1; i <= cyield->GetXaxis()->GetNbins(); i++ ){
			float v = cyield->GetBinContent( i );
			float e = cyield->GetBinError( i );
			float pt = cyield->GetBinCenter( i );

			float vTBW = hTBW->GetBinContent( hTBW->GetXaxis()->FindBin(pt) );

			yratio->SetBinContent( i, v / vTBW );
			yratio->SetBinError( i, e / vTBW );
		}
		
		gPad->SetLogy(0);
		rpl.style( yratio ).set( config, "style.yratio" ).draw();
		rp.next();


		rpl.style( "mass" ).set( "style.massfit" ).draw();
		rp.next();
		rpl.style( "width" ).set( "style.width" ).draw();
		rp.next();
		rpl.style( "SoverB" ).set( "style.SoverB" ).draw();
		rp.next();
		rpl.style( "significance" ).set( "style.significance" ).draw();
		rp.next();
	}

};



#endif
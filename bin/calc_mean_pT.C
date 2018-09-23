


double calc_mean_pT( double low, double high) {

	TFile * f = new TFile( "/Users/jdb/bnl/work/dimuonAna/data/Cocktail/MTDMTD-phi_mumu.root" );

	TF1 * f1 = (TF1*)f->Get( "phi_mumu/phi_pT" );

	double wm = f1->Mean( low, high );
	double m = (low+high) / 2.0;

	cout << "meanPt = " << wm << " unweighted = " << m << endl;

	return wm;
}
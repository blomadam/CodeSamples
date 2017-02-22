#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

class DRTableBP {
	private: 
		double iGridLowerBound , iGridUpperBound , iGridStep ; // k_in grid definitions
		double oGridLowerBound , oGridUpperBound , oGridStep ; // k_out grid definitions
		double aGridLowerBound , aGridUpperBound , aGridStep ; // k_th grid definitions
		double tGridLowerBound , tGridUpperBound , tGridStep ; // theta grid definitions
		double pGridLowerBound , pGridUpperBound , pGridStep ; // phi grid definitions
		int Ignum, Ognum, Agnum, Tgnum, Pgnum; // number of grid points in each variable
		double *Ipts ,*Opts ,*Apts , *Tpts, *Ppts ;  // arrays to hold the grid
		double *gsig5;   // xsecs from table
	public: 
		void FillFromFile(); // reads data file and fills arrays listed below
		double CalcSigma(double k_in,double k_out, double k_th, 
				 double theta, double phi ); // takes kinematics and returns 5F CM xsecs 
		double cubicInterpolate (double [4], double ) ;
		double nCubicInterpolate (int n, double* p, double coordinates[]);
};

void DRTableBP::FillFromFile() {
	std::cout << "reading xsec file\n";
	std::ifstream input("VCSlookup.dat");
	
	double  a,b,c,d,e,f,g,h,i,m,n,o,p,q,r;
	
//	read first line of lookup table for structure information
// the parameter bounds and step size are listed and used to great the arrays to hold
// the cross section grid in each variable.
        input >> a  >> b  >> c  >> d  >> e  >> f  >> g  >> h >> i >> m >>n >>o >>p >>q >>r ;

		iGridLowerBound = a;
		iGridUpperBound = b;
		iGridStep = c;
		oGridLowerBound = d;
		oGridUpperBound = e;
		oGridStep = f;
		aGridLowerBound = g;
		aGridUpperBound = h;
		aGridStep = i;
		tGridLowerBound = m;
		tGridUpperBound = n;
		tGridStep = o;
		pGridLowerBound = p;
		pGridUpperBound = q;
		pGridStep = r;

        Ignum = (int)((iGridUpperBound - iGridLowerBound) / iGridStep +0.5);
        Ognum = (int)((oGridUpperBound - oGridLowerBound) / oGridStep +0.5);
        Agnum = (int)((aGridUpperBound - aGridLowerBound) / aGridStep +0.5);
        Tgnum = (int)((tGridUpperBound - tGridLowerBound) / tGridStep +0.5);
        Pgnum = (int)((pGridUpperBound - pGridLowerBound) / pGridStep +0.5);

		Ipts = new double[ Ignum+1 ];
		Opts = new double[ Ognum+1 ];
		Apts = new double[ Agnum+1 ];
		Tpts = new double[ Tgnum+1 ];
		Ppts = new double[ Pgnum+1 ];
		gsig5 = new double[(Ignum+1) * (Ognum+1) * (Agnum+1) * (Tgnum+1) * (Pgnum+1) ];

	int Ipos,Opos,Apos, Tpos, Ppos, sigpos; 
	while (input >> a  >> b  >> c   >> e  >> f >> g  )
	{
// Fix positions for order in the table
		Ipos = (int)((a-iGridLowerBound)/iGridStep+0.5);
		Opos = (int)((b-oGridLowerBound)/oGridStep+0.5);
		Apos = (int)((c-aGridLowerBound)/aGridStep+0.5);
		Tpos = (int)((e-tGridLowerBound)/tGridStep+0.5);
		Ppos = (int)((f-pGridLowerBound)/pGridStep+0.5);
// Fix sigpos for the new layout
		sigpos = Ipos * (Ognum+1) * (Agnum+1) * (Tgnum+1) * (Pgnum+1)
				+ Opos * (Agnum+1) * (Tgnum+1) * (Pgnum+1)
				+ Apos * (Tgnum+1) * (Pgnum+1)
				+ Tpos * (Pgnum+1)
				+ Ppos  ;
// store the parameter values and cross section in corresponding array
		Ipts[Ipos] = a;
		Opts[Opos] = b;
		Apts[Apos] = c;
		Tpts[Tpos] = e;
		Ppts[Ppos] = f;
		gsig5[ sigpos ] = g;
// cout << sigpos << " " << g << endl;
		
	}
		std::cout << "finished reading xsec file\n";
	return;
}



// calculate cross section
double DRTableBP::CalcSigma(
		double k_in, // GeV
		double k_out, // GeV
		double k_th, // degrees
		double theta, // degrees
		double phi  // degrees
		) {
	
	//double const R2D = 180.0/3.14159265;      // convert radians to degrees
	int Ival = int ( floor ( (k_in  - iGridLowerBound) / iGridStep  )  );
	int Oval = int ( floor ( (k_out - oGridLowerBound) / oGridStep  )  );
	int Aval = int ( floor ( (k_th - aGridLowerBound) / aGridStep  )  );
	int Tval = int ( floor ( (theta - tGridLowerBound) / tGridStep  )  );
	int Pval = int ( floor ( (phi   - pGridLowerBound) / pGridStep  )  );
	
	double p[4][4][4][4][4];  // 6 dimensional lookup
	double sigma5;
		
	
	if (Ival < 0 || Oval < 0 || Aval < 0 || Tval < 0  || Pval < 0 
	     || Ival > Ignum  || Oval > Ognum  || Aval > Agnum  || Tval > Tgnum || Pval > Pgnum 

	     ) {
/*		cout << k_in << " " << "Ival: " << Ival << " of " << Ignum << endl;
        cout << k_out << " "<< "Oval: " << Oval << " of " << Ognum << endl;
        cout << k_th << " " << "Aval: " << Aval << " of " << Agnum << endl;
        cout << theta << " "<< "Tval: " << Tval << " of " << Tgnum << endl;
        cout << phi << " "  << "Pval: " << Pval << " of " << Pgnum << endl;
		cout << "look up out of table bounds\n";
*/
			return 0.0 ;
	}
	
	
	for (int ii = -1; ii < 3; ii++) {
	  for (int oo = -1; oo < 3; oo++) {
	    for (int aa = -1; aa < 3; aa++) {
	      for (int tt = -1; tt < 3; tt++) {
	      for (int pp = -1; pp < 3; pp++) {
		 
		 // code to mirror parameter values around boundary limits
		 // need to include grid previous parameters Tpos = (Tval+tt) *(Pgnum+1)
		 // Then swap into testval below
//FIXME: note the issues for inner bound if tables do not run full width 0-180 deg
		 int Ppos = Pval+pp ;  
		 if (Pval+pp < 0) {
		 Ppos = 0-(Pval+pp);
//		 cout << Ppos  << " " << pp<< endl;
		 }
		 else if (Pval+pp > Pgnum) {
		  Ppos = 2*Pgnum - (Pval+pp);
//		 cout << Ppos  << " " << pp<< endl;
		  }
			  
			  int testval = (Ival+ii) *(Ognum+1)*(Agnum+1)*(Tgnum+1)*(Pgnum+1) 
				+ (Oval+oo) *(Agnum+1)*(Tgnum+1)*(Pgnum+1) 
				+ (Aval+aa) *(Tgnum+1)*(Pgnum+1) 
				+ (Tval+tt) *(Pgnum+1) 
				+ Ppos ;
//				cout << testval  << " " << gsig5[testval]<< endl;
			  if ( testval<0 || testval  > (Ignum+1) * (Ognum+1) * (Agnum+1) * (Tgnum+1) * (Pgnum+1) ){
//				cout << testval <<endl;
//				cout << "sig lookup out of bounds\n";
			return 0.0 ;
				}
			  
// Fix for proper order
		 p[ii+1][oo+1][aa+1][tt+1][pp+1] = 
			gsig5[testval] ;
	      }  // close pp			
	      }  // close tt			
	    } // close aa
	  } // close oo
	} // close ii
				  
	double pts[5] = {(k_in  - Ipts[Ival]) / iGridStep,
			(k_out  - Opts[Oval]) / oGridStep,
			(k_th   - Apts[Aval]) / aGridStep,
			(theta  - Tpts[Tval]) / tGridStep, 
			(phi    - Ppts[Pval]) / pGridStep };
	sigma5=nCubicInterpolate(5,(double*)p,pts);
	
//	cout << sigma5 <<endl;
	

	return sigma5;

}




/*------------------------------------------------------------------------------------
  The following routines are adapted from the algorithms
  presented at http://www.paulinternet.nl/?page=bicubic
------------------------------------------------------------------------------------*/
double DRTableBP::cubicInterpolate (double p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2]
	        - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}
double DRTableBP::nCubicInterpolate (int n, double* p, double coordinates[]) {
	assert(n > 0);
	if (n == 1) {
		return cubicInterpolate(p, *coordinates);
	}
	else {
		double arr[4];
		int skip = 1 << (n - 1) * 2;
		arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1);
		arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1);
		arr[2] = nCubicInterpolate(n - 1, p + 2*skip, coordinates + 1);
		arr[3] = nCubicInterpolate(n - 1, p + 3*skip, coordinates + 1);
		return cubicInterpolate(arr, *coordinates);
	}
}

int main () {

// define a lookup table instance
  class DRTableBP *BPtablexsecs;
	BPtablexsecs = new DRTableBP;
// fill arrays with data from file
	  BPtablexsecs->FillFromFile();
// run a test point... should get 497
	std::cout << BPtablexsecs->CalcSigma(1.096, 0.644,31.1,130.5,166.2) << std::endl;
	std::cout << BPtablexsecs->CalcSigma(1.096, 0.654,30.1,138.3,176.8) << std::endl;


}


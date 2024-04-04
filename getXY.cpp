// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "index.h"

// these are for lattice matched material

double Band(double x, double y)
{
	double Eg=0.36+2.093*y+0.629*x+0.577*y*y+0.436*x*x+1.013*x*y-2.0*x*y*(1-x-y);

	return (Eg);
}

double MondryBand(double x)
{
	double Eg=0.75+0.72*x;
	return(Eg);
}

double GetXmax(void)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP;

	double xmax=(a-a_InAs)/(a_GaAs-a_InAs );

	return(xmax);
}

double GetYmax(void)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP;

	double ymax=(a-a_InAs)/(a_AlAs-a_InAs );

	return(ymax);
}

double GetX(double y)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP;

	double x=(a-a_AlAs*y-a_InAs*(1-y))/(a_GaAs-a_InAs );

	return(x);
}

double GetY(double x)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP;

	double xmax=(a-a_InAs)/(a_GaAs-a_InAs );

	if (x>xmax) return(0);

	double y=(a-a_GaAs*x-a_InAs*(1-x))/(a_AlAs-a_InAs);

	return(y);
}

double GetX(double Strain, double y)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP/(1+Strain);

	double x=(a-a_AlAs*y-a_InAs*(1-y))/(a_GaAs-a_InAs );

	return(x);
}

double GetXmax(double Strain)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP/(1+Strain);

	double xmax=(a-a_InAs)/(a_GaAs-a_InAs );

	return(xmax);
}

double GetYmax(double Strain)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP/(1+Strain);

	double ymax=(a-a_InAs)/(a_AlAs-a_InAs );

	return(ymax);
}

double GetY(double Strain, double x)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A
	double a=a_InP/(1+Strain);

	double xmax=(a-a_InAs)/(a_GaAs-a_InAs );

	if (x>xmax) return(0);

	double y=(a-a_GaAs*x-a_InAs*(1-x))/(a_AlAs-a_InAs);

	return(y);
}

double Strain(double x, double y)
{
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A

	//Lattice constant for InGaAlAs
	double a=a_InAs*(1-x-y)+a_GaAs*x+a_AlAs*y; //A

 	// Strain
	double epsilon=(a_InP-a)/a;

	return(epsilon);
}

void GetXY(double SetBand, double SetStrain, double &x, double &y)
{
	double xmax=GetXmax();
	double xmin=0.0;
	double B,Bmin,Bmax;

	// check for out of range
	if (SetBand>(Bmax=Band(x=xmin,y=GetY(SetStrain,xmin)))) {
		return;
	}
	else if (SetBand<(Bmin=Band(x=xmax,y=GetY(SetStrain,xmax)))) {
		return;
	}

	do {
		x=(Bmax-SetBand)/(Bmax-Bmin)*(xmax-xmin)+xmin;
		y=GetY(SetStrain,x);
		if ((B=Band(x,y))>SetBand) {
			Bmax=B;
			xmin=x;
		}
		else {
			Bmin=B;
			xmax=x;
		}
	}
	while (fabs(B-SetBand)>1E-10);

}

void GetXY(double SetBand, double &x)
{
	x=(SetBand-0.75)/0.72;
}

void GetXY(double SetBand, double &x, double &y)
{
	double xmax=GetXmax();
	double xmin=0.0;
	double B,Bmin,Bmax;

	// check for out of range
	if (SetBand>(Bmax=Band(x=xmin,y=GetY(xmin)))) {
		return;
	}
	else if (SetBand<(Bmin=Band(x=xmax,y=GetY(xmax)))) {
		return;
	}

	do {
		x=(Bmax-SetBand)/(Bmax-Bmin)*(xmax-xmin)+xmin;
		y=GetY(x);
		if ((B=Band(x,y))>SetBand) {
			Bmax=B;
			xmin=x;
		}
		else {
			Bmin=B;
			xmax=x;
		}
	}
	while (fabs(B-SetBand)>1E-10);

}

void Material(double x, double y) {

	cout <<"x=" << x << "\ty=" << y << endl;
	// Parameters
	double Eg_GaAs=1.424; //eV
	double Eg_AlAs=2.168; //eV
	double Eg_InAs=0.354; //eV


	// Band gap
	double Eg=0.36+2.093*y+0.629*x+0.577*y*y+0.436*x*x+1.013*x*y-2.0*x*y*(1-x-y);
	cout <<"bandgap" << endl;
	cout <<"Eg=" <<Eg<< endl;

	// Band gap with strain
	// Lattice constant;
	double a_GaAs=5.6533; //A
	double a_InAs=6.0584; //A
	double a_AlAs=5.660; //A

	double a_InP=5.8688; //A

	//Hydrostatic deformation potential for conduction band
	double ac_GaAs=-7.17;  //(eV)
	double ac_InAs=-5.08;  //(eV)
	double ac_AlAs=-5.64;  //(eV)

	double ac_InP=-5.04; //(eV)
	//Hydrostatic deformation potential for valence band
	double av_GaAs=1.16;   //(eV)
	double av_InAs=1.00;   //(eV)
	double av_AlAs=2.47;   //(eV)

	double av_InP=1.27;   //(eV)
	//Shear deformation potential for valence band
	double b_GaAs=-1.7; //(eV)
	double b_InAs=-1.8; //(eV)
	double b_AlAs=-1.5; //(eV)

	double b_InP=-1.7; //(eV)
	//Elastic Stiffness constants
	double C11_GaAs=11.879; //10^11 dyn/cm^2
	double C12_GaAs=5.376; //10^11 dyn/cm^2

	double C11_InAs=8.329; //10^11 dyn/cm^2
	double C12_InAs=4.526; //10^11 dyn/cm^2

	double C11_AlAs=12.5;  //10^11 dyn/cm^2
	double C12_AlAs=5.34;  //10^11 dyn/cm^2

	double C11_InP=10.11; //10^11 dyn/cm^2
	double C12_InP=5.61;  //10^11 dyn/cm^2

	//Lattice constant for InGaAlAs
	double a=a_InAs*(1-x-y)+a_GaAs*x+a_AlAs*y; //A

	//Hydrostatic deformation potential for conduction band for InGaAlAs
	double ac=ac_InAs*(1-x-y)+ac_GaAs*x+ac_AlAs*y; //(eV)

	//Hydrostatic deformation potential for valence band for InGaAlAs
	double av=av_InAs*(1-x-y)+av_GaAs*x+av_AlAs*y; //(eV)

	//Shear deformation potential for valence band for InGaAlAs
	double b=b_InAs*(1-x-y)+b_GaAs*x+b_AlAs*y; //(eV)

	//Elastic Stiffness constants for InGaAlAs
	double C11=C11_InAs*(1-x-y)+C11_GaAs*x+C11_AlAs*y; //10^11 dyn/cm^2
	double C12=C12_InAs*(1-x-y)+C12_GaAs*x+C12_AlAs*y; //10^11 dyn/cm^2

	// Effective mass
	//Electron effective mass:
	double me_GaAs=0.067; //(m0)
	double me_InAs=0.023; //(m0)
	double me_AlAs=0.15; //(m0)

	double me=me_InAs*(1-x-y)+me_GaAs*x+me_AlAs*y;

	//HH effective mass
	double mhh_GaAs=0.50; //(m0)
	double mhh_InAs=0.40; //(m0)
	double mhh_AlAs=0.79; //(m0)

	double mhh=mhh_InAs*(1-x-y)+mhh_GaAs*x+mhh_AlAs*y;

	cout << "The effective mass" << endl;
	cout << " for electrons       for holes" << endl;
	cout << "m_e/m_0=" <<me<<";           m_hh/m_0=" <<mhh<<endl;

 	// Strain
	double epsilon=(a_InP-a)/a_InP;
	double epsilon_zz=-2*C12/C11*epsilon;

	cout <<"\nThe strain in the plane of the epitaxial growth is:" << endl;
	cout <<"epsilon="<< epsilon << endl;

	cout << "\nThe strain in the perpendicular direction is:" << endl;
	cout << "epsilon_zz=" << epsilon_zz << endl;
	
	double P_epsilon=-2*av*(1-C12/C11)*epsilon;
	double Q_epsilon=-b*(1+2*C12/C11)*epsilon;

	//double deltaEc=2*ac*(1-C12/C11)*epsilon;
	//double deltaEhh=-P_epsilon-Q_epsilon;
	//double deltaElh=-P_epsilon+Q_epsilon;
	double deltaEc=2*ac*(1-C12/C11)*epsilon+P_epsilon;
	double deltaEhh=-Q_epsilon;
	double deltaElh=Q_epsilon;

	cout << "\nThe shift in the conduction band" << endl;
	cout << "delta Ec=" << deltaEc <<";       delta Ehh=" << deltaEhh << ";        delta Elh=" << deltaElh << endl;

	// The strained bandgap;
	double Ec_hh=Eg+deltaEc-deltaEhh;
	double Ec_lh=Eg+deltaEc-deltaElh;

	cout << "\nThe strained bandgap" << endl;
	cout << "Ec_hh=" << Ec_hh << ";        Ec_lh=" << Ec_lh << endl;

	// Band Offset
	// Model-Solid Theory

	//average valence band position
	double Evav_GaAs=-6.92; //eV
	double Evav_InAs=-6.67; //eV
	double Evav_AlAs=-7.49; //eV

	double Evav_InP=-7.04;  //(eV)

	//Spin orbit split-off energy
	double Delta_GaAs=0.34; //eV
	double Delta_InAs=0.38;  //eV
	double Delta_AlAs=0.28;  //eV

	double Delta_InP=0.11; //(eV)

	// the average valence subband energy
	double Evav=Evav_InAs*(1-x-y)+Evav_GaAs*x+Evav_AlAs*y;
	// the spin-orbit split-off band energy
	double Delta=Delta_InAs*(1-x-y)+Delta_GaAs*x+Delta_AlAs*y;

	// The valence band position of a quaternary
	// for hh (compresive strain)
	double Ev_hh=Evav+Delta/3+deltaEhh;
	// for lh (tensile strain)
	double Ev_lh=Evav+Delta/3+deltaElh;

	cout << "\nThe valence bands position in Model-Solid Theory:" << endl;
	cout << "compresive strain;    tensile strain" << endl;
	cout << "Ev_hh=" << Ev_hh << ";           Ev_lh=" << Ev_lh << endl;

	// The conduction band position may be calculated by simply adding
	// the strained bandgap energy to the valence band position.

	double Echh=Ev_hh+Ec_hh; // for hh (compresive strain)
	double Eclh=Ev_lh+Ec_lh; // for lh (tensile strain)


	cout << "\nThe conduction band position in Model-Solid Theory:" << endl;
	cout << " compresive strain;    tensile strain" << endl;
	cout << "Echh=" << Echh << ";           Eclh=" << Eclh << endl;


	// InP parameters
	//InP band gap
	double Eg_InP=1.344; //(eV)

	double P_epsilon_InP=-2*av_InP*(1-C12_InP/C11_InP)*epsilon;
	double Q_epsilon_InP=-b_InP*(1+2*C12_InP/C11_InP)*epsilon;

	double deltaEc_InP=2*ac_InP*(1-C12_InP/C11_InP)*epsilon;
	double deltaEhh_InP=-P_epsilon_InP-Q_epsilon_InP;
	double deltaElh_InP=-P_epsilon_InP+Q_epsilon_InP;

	// The strained bandgap for InP;
	double Ec_hh_InP=Eg_InP+deltaEc_InP-deltaEhh_InP;
	double Ec_lh_InP=Eg_InP+deltaEc_InP-deltaElh_InP;

	// The valence band position of a InP
	//for hh (compresive strain)
	double Ev_hh_InP=Evav_InP+Delta_InP/3+deltaEhh_InP;
	// for lh (tensile strain)
	double Ev_lh_InP=Evav_InP+Delta_InP/3+deltaElh_InP;


	// hh
	double Evw=Ev_hh; // the valence band positions in the well
	double Evb=Ev_hh_InP; //  the valence band positions in barrier materials,
	double Egw=Ec_hh; // strain adjusted band gaps for the well
	double Egb=Ec_hh_InP; // strain adjusted band gaps for the barrier.

	double DeltaEc_DeltaEg_hh=1-(Evw-Evb)/(Egb-Egw);
	//lh
	double Evw_lh=Ev_lh; // the valence band positions in the well
	double Evb_lh=Ev_lh_InP; //  the valence band positions in barrier materials,
	double Egw_lh=Ec_lh; // strain adjusted band gaps for the well
	double Egb_lh=Ec_lh_InP; // strain adjusted band gaps for the barrier.

	double DeltaEc_DeltaEg_lh=1-(Evw_lh-Evb_lh)/(Egb_lh-Egw_lh);

	cout << "\nThe conduction band-edge discontinuity in Model-Solid Theory:" << endl;
	cout << " compresive strain;    tensile strain" << endl;
	cout << "(Delta Ec)/(Delta Eg)_hh=" << DeltaEc_DeltaEg_hh << ";           (Delta Ec)/(Delta Eg)_lh=" << DeltaEc_DeltaEg_lh << endl;

	// Band Offset
	// Harrison’s Model

	// conduction band position
	double EcH_GaAs=1.53;  //eV
	double EcH_InAs=0.801; //eV
	double EcH_AlAs=2.5255;  //eV

	double EcH_InP=1.35; //(eV)

	//valence band position
	double EvH_GaAs=0.111;  //eV
	double EvH_InAs=0.441; //eV
	double EvH_AlAs=-0.4245;  //eV

	double EvH_InP=0;

	// conduction band position of InGaAlAs
	double EcH=EcH_InAs*(1-x-y)+EcH_GaAs*x+EcH_AlAs*y;
	//valence band position f InGaAlAs
	double EvH=EvH_InAs*(1-x-y)+EvH_GaAs*x+EvH_AlAs*y;

	double Ev_H_hh=EvH+deltaEhh; //for hh (compresive strain)
	double Ev_H_lh=EvH+deltaElh; //for lh (tensile strain)

	double Ec_H= EcH+deltaEc;

	cout << "\nThe conduction and valece bands in Harrison Model:" << endl;
	cout << "Concuction band;      VB-compresive strain;    VB-tensile strain" << endl;
	cout << "Ec_H=" << Ec_H << ";           Ev_H_hh=" << Ev_H_hh << ";           Ev_H_lh=" << Ev_H_lh << endl;


	// The conduction band-edge discontinuity
	double Ev_H_hh_InP=EvH_InP+deltaEhh_InP; //for hh (compresive strain)
	double Ev_H_lh_InP=EvH_InP+deltaElh_InP; //for lh (tensile strain)

	double Ec_H_InP= EcH_InP+deltaEc_InP;

	//%%%%%%%%
	double EcH_b=Ec_H_InP;
	double EvH_b_hh=Ev_H_hh_InP;
	double EvH_b_lh=Ev_H_lh_InP;

	double EcH_w=Ec_H;
	double EvH_w_hh=Ev_H_hh;
	double EvH_w_lh=Ev_H_lh;
	double DeltaEc_DeltaEg_H_hh=(EcH_b-EcH_w)/((EvH_w_hh-EvH_b_hh)+(EcH_b-EcH_w));
	double DeltaEc_DeltaEg_H_lh=(EcH_b-EcH_w)/((EvH_w_lh-EvH_b_lh)+(EcH_b-EcH_w));

	cout << "\nThe conduction band-edge discontinuity in Harrison Model" << endl;
	cout << " compresive strain;    tensile strain" << endl;
	cout << "(Delta Ec)/(Delta Eg)_hh=" << DeltaEc_DeltaEg_H_hh << ";           (Delta Ec)/(Delta Eg)_lh=" << DeltaEc_DeltaEg_H_lh<< endl;

}

int _tmain(int argc, _TCHAR* argv[])
{
	double x,y,z;
	double BG;
	double Eg,Str;

/*	//Material(x,y);
	GetXY(0.8,-0.005,x,y);
	cout << "Eg(" << x << "," << y << ") = " << Band(x,y) << "\teps=" << Strain(x,y) << endl;
	GetXY(1.0,0.0,x,y);
	cout << "Eg(" << x << "," << y << ") = " << Band(x,y) << "\teps=" << Strain(x,y) << endl;
	GetXY(1.2,0.005,x,y);
	cout << "Eg(" << x << "," << y << ") = " << Band(x,y) << "\teps=" << Strain(x,y) << endl;

	x=0.1831;  y=0.283;
	cout << "Eg(" << x << "," << y << ") = " << Band(x,y) << "\teps=" << Strain(x,y) << endl;
	x=0.2409;  y=0.187;
	cout << "Eg(" << x << "," << y << ") = " << Band(x,y) << "\teps=" << Strain(x,y) << endl;
	x=0.2412;  y=0.29;
	cout << "Eg(" << x << "," << y << ") = " << Band(x,y) << "\teps=" << Strain(x,y) << endl;

	fstream fs;
	fs.open ("test.txt", fstream::in | fstream::out | fstream::trunc);

	fs << "x\tEg\tEg\teps" << endl;
	for (int i=0;i<101;i++) {
		z=0.01*(double)i;
		x=0.47*(1-z);
		y=0.48*z;
		fs << z << "\t" << MondryBand(z) << "\t" << Band(x,y) << "\t" << Strain(x,y) << endl;
	}

	fs.close();*/
	cin >> BG;
	GetXY(BG,0.000,x,y);
	cout << "x=" << x << "\ty=" << y <<endl;;
	
	//~ x=0.2409;  y=0.187;
	//~ Material(x,y);

	// for InGaAs
	//Material(1.0,0.0);

	return 0;
}


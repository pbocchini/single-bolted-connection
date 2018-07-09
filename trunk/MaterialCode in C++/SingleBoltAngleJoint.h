/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written by Liyang Ma (lim215@lehigh.edu)
// Created 2018-02-01 $
// 
//Description: This file contains the class definition for SingleBoltAngleJoint 
//SingleBoltAngleJoint describes the hysterestic beahvior of single bolted angle joint, it needs to be applied to zero-length element
                                                                        
#ifndef SingleBoltAngleJoint_h
#define SingleBoltAngleJoint_h



#include <UniaxialMaterial.h>
#include <vector>
#include <algorithm>


typedef std::vector<int>    int_vec;
typedef std::vector<double> double_vec;

class SingleBoltAngleJoint : public UniaxialMaterial
{
public:
	SingleBoltAngleJoint( int tag,
		//Slip load
		double Pslip,
		//Plate geometry properties
		double width_brace, double thickness_brace, double thickness_leg, double  Bolt_diameter, double Hole_diameter, double e1,
		//Plate material properties
		double Fy_brace, double Fu_brace, double Fy_leg,
		double E, double v
		);
	SingleBoltAngleJoint(void);
	virtual ~SingleBoltAngleJoint();

	const char *getClassType(void) const { return "SingleBoltAngleJoint"; };


	//functions
	double Line(double x, double k, double u, double f);
	double Intersectu(double k1, double u1, double f1, double k2, double u2, double f2);
	double Intersectf(double k1, double u1, double f1, double k2, double u2, double f2);
	double Bisectionu(double u, double f, double Kb, double_vec REP, double u_b0, double R, double Ki, int status, double f_b0, double Kslip);
	double Bisectionf(double u, double f, double Kb, double bisectionu);
	double RichardEquation(double_vec REP, double f_b0, double Kslip, double R, double Ki, double u, double u_b0);

	void CalcBreakpoints(int dir, int status, double Pslip, double u_bc, double f_bc, double u_bt, double f_bt, double u_b0c, double f_b0c,
		double u_b0t, double f_b0t, double u_s0c, double f_s0c, double u_s0t, double f_s0t, double u_current, double f_current,
		double Keb, double Kti, double Kci, double_vec REPc, double_vec REPt, double Rc, double Rt, double Kslip, double u_transition, double f_transition,
		double_vec& ubkp, double_vec& fbkp, double_vec& kbkp);
	void CalcForceStatusDir(double_vec ubkp, double_vec fbkp, double_vec kbkp, double u_trial, int status, int dir, double u_current, double f_current, double u_st, double f_st,
		double u_sc, double f_sc, double Kslip, double Kti, double Kci, double u_transition, double f_transition,
		double u_bc, double f_bc, double u_bt, double f_bt, double_vec REPt, double_vec REPc, double Rt, double Rc, double u_b0t, double f_b0t, double u_b0c, double f_b0c, double Keb, double Pslip,
		double& f_trial, int& status_trial, int& dir_trial, double& u_transition_trial, double& f_transition_trial, double& u_bc_trial, double& f_bc_trial, double& u_bt_trial,
		double& f_bt_trial, double& u_b0t_trial, double& f_b0t_trial, double& u_b0c_trial, double& f_b0c_trial, double& k_trial);
	// methods must implement for the class to link successfully with OpenSees
	UniaxialMaterial *getCopy(void);
	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	// MATERIAL INPUTS ----------------------------------------------------------------------------
	double Pslip; // Slip load
	double width_brace; // width of the brace plate
	double thickness_brace; //thickness of the brace plate
	double thickness_leg; // thickness of the leg plate 
	double Bolt_diameter; // diameter of the bolt
	double Hole_diameter; // diameter of the bolt hole
	double e1; // distance from the edge of the plate to the center of the bolt hole
	double Fy_brace; //yield strength of the brace plate
	double Fu_brace; //ultimate strength of the brace plate
	double Fy_leg;   //yield strength of the leg plate
	double E; // Young's modulus
	double v; // Poisson ratio

   // FIXED PROPERTIES -----------------------------------------------------------------------------
	double clearance; // difference between bolt diameter and hole diameter
	double G; //shear modulus 
	double Leff; //effective length calulated using euro code 3 with 30 degrees force angle
	double An; //area of the cross section
	double Aeff; //half area of the cross section
	double Af; //area calculated as width times thickness
	double Anet; //area calculated as (width-bolt diameter)*thickness
	double Keb; //Elastic stiffness of the brace plate using Leff caculated as Aeff*E/Leff
   // INITIAL COMPRESSIVE STIFFNESS
	double Kbr; //brace plate bearing stiffness
	double Knet; //brace plate net section stiffness
	double Abolt; //area of bolt
	double Lbolt; // length of the bolt
	double Ibolt; //moment of inertia of the bolt
	double Kboltshear; //bolt shear stiffness
	double Kboltbending; //bolt bending stiffness
	double KbrLeg; //leg plate bearing stiffness
	double Kci; //initial compressive stiffness
   // INITIAL TENSILE STIFFNESS
	double Kbben; //brace plate bending stiffness
	double Kbv; //brace plate shear stiffness
	double Kti; //initial tensile stiffness
   // COMPRESSIVE CAPACITY 
	double Rc;
   // TENSILE CAPACITY
	double Rt;
   // RICHARD EQUATION PARAMETER
	double_vec REPc;   //RE parameter for compression
	double_vec REPt;   //RE parameter for tension
   // OTHERS
	double Kslip; // Assign very small stiffness to slippage to avoid numerical problems
	double u_inc; // u increment used for tangent evaluation 


   
   // VARIBLES SUBJECT TO CHANGE-------------------------------------------------------------------
   // INITIAL COMPRESSIVE BEARING DISPLACEMENT 
	double u_s0c;  // compression slipping displacment 
	double f_s0c;

	double u_b0c;  // start of Richard equation
	double f_b0c; 
	double u_b0c_trial; 
	double f_b0c_trial;

   // INITIAL TENSILE BEARING DISPLACEMENT 
	double u_s0t;  // tensile slipping displacment 
	double f_s0t;

	double u_b0t;  // start of Richard equation
	double f_b0t;
	double u_b0t_trial; 
	double f_b0t_trial;

   // bearing displacment and force
	double u_bc; // start of bearing; different from start of Richard Equation 
	double f_bc;
	double u_bc_trial;
	double f_bc_trial;

	double u_bt;
	double f_bt;
	double u_bt_trial; 
	double f_bt_trial;
   // slipping displacement and force 
	double u_sc;
	double f_sc;
	double u_st;
	double f_st;
 
  // break points vector
	double_vec ubkp;  //6*1 vector for u_current<=ubkp[1] and u_current>=ubkp[6] we use Richard Equation 
	double_vec fbkp;  //6*1 vector
	double_vec kbkp;  //5*1 vector

  // status 1 tension bearing nonstick -1 compression bearing nonstick 
  // status 4 tension bearing stick - 4 compressiong bearing stick
  // status 0 before bearing
	int status;
	int status_trial;
  
  //direction 1 loading - 1 reversal loading 0 for start
	int dir;
	int dir_trial;
  
  // CURRENT COMMITTED U  F  K
	double u_current;
	double f_current;
	double k_current;
  // TRIAL U F K
	double u_trial; // can be defined as trialStrain from setTrailStrain
	double f_trial;
	double k_trial;
  // TRANSITION U F FOR BEARING STICK PHASE
	double u_transition;
	double f_transition;
	double u_transition_trial;
	double f_transition_trial;

};
#endif




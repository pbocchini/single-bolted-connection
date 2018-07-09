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

#define _USE_MATH_DEFINES   //never seen beingi used in OPENSEES
#include <stdlib.h>
#include <SingleBoltAngleJoint.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <elementAPI.h>
#include <math.h>
#include <algorithm>

void*
OPS_SingleBoltAngleJoint(void)
{   
	//Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;
	//
	//Parse the input line for the material parameters
	//
	int    iData[1];
	double dData[12];
	int numData;
	numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SingleBoltAngleJoint tag" << endln;
		return 0;
	}

	numData = 12;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING not 12 input parameters \n";
		return 0;
	}
	//                                                                                                                                           
	// create a new material                                                                                                                     
	//                                                                                                                                           
	//Allocate the material
	theMaterial = new SingleBoltAngleJoint(iData[0], 
		                                   dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		                                   dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);
	// in case there was a problem with material creation
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SingleBoltAngleJoint\n";
		return 0;
	}
	return theMaterial;
}

SingleBoltAngleJoint::SingleBoltAngleJoint(int tag,
	//Slip load
	double _Pslip,
	//Plate geometry properties
	double _width_brace, double _thickness_brace, double _thickness_leg, double  _Bolt_diameter, double _Hole_diameter, double _e1,
	//Plate material properties
	double _Fy_brace, double _Fu_brace, double _Fy_leg,
	double _E, double _v
) :
	UniaxialMaterial(tag, 0), 
	Pslip(_Pslip), 
	width_brace(_width_brace), thickness_brace(_thickness_brace), thickness_leg(_thickness_leg), Bolt_diameter(_Bolt_diameter), Hole_diameter(_Hole_diameter), e1(_e1),
	Fy_brace(_Fy_brace), Fu_brace(_Fu_brace), Fy_leg(_Fy_leg),
	E(_E), v(_v)
	
{
	clearance = Hole_diameter - Bolt_diameter;
	G = E / (2 * (1 + v));
	//Effective length calulated using euro code 3 with 30 degrees force angle
	Leff = pow((width_brace + (width_brace - thickness_brace) / 2) * 3, 0.5);
	//Elastic stiffness of the brace plate using Leff
	An = width_brace*thickness_brace + (width_brace - thickness_brace)*thickness_brace;
	Aeff = An / 2;
	Af = width_brace*thickness_brace;
	Anet = (width_brace - Bolt_diameter)*thickness_brace;
	Keb = Aeff*E / Leff;
	//Initial Compressive stiffness
	Kbr = 120 * Fy_brace*thickness_brace*Bolt_diameter;  //braceplate bearing stiffness
	Knet = Keb;
	Abolt = pow(M_PI*(Bolt_diameter / 2), 2);
	Lbolt = thickness_brace + thickness_leg;
	Ibolt = pow(M_PI / 4 * (Bolt_diameter / 2), 4);
	Kboltshear = Abolt*G / Lbolt * 32 / 37;
	Kboltbending = 3 * E*Ibolt / (pow(Lbolt , 3));
	KbrLeg = 120 * Fy_leg*thickness_leg*Bolt_diameter;
	Kci = 1 / (1 / Kbr + 1 / Knet + 1 / Kboltshear + 1 / Kboltbending + 1 / KbrLeg);
	//Initial Tensile stiffness
	Kbben = 32 * E*thickness_brace*pow((e1 / Bolt_diameter - 1 / 2) , 3); //brace plate bending stiffness
	Kbv = 6.67*G*thickness_brace*(e1 / Bolt_diameter - 1 / 2); //brace plate shear stiffness
	Kti = 1 / (1 / Kci + 1 / Kbben + 1 / Kbv);
	//Compressive Capacity   
	Rc = ((2.4*Bolt_diameter*thickness_brace*Fu_brace) < (Af*Fy_brace)) ? 2.4*Bolt_diameter*thickness_brace*Fu_brace : Af*Fy_brace;
	//Tensile Capacity
	Rt = (Anet*Fu_brace < 0.7*Fu_brace * 2 * thickness_brace*(e1 - Hole_diameter / 2) / cos(M_PI / 6)) ? Anet*Fu_brace : 0.7*Fu_brace * 2 * thickness_brace*(e1 - Hole_diameter / 2) / cos(M_PI / 6);
	//Richard Equation Parameter :
	REPc[0] = 7.28899865552825;   //Compression
	REPc[1] = 0.00711376636993;
	REPc[2] = 2.77499874147550;
	REPc[3] = 0.33019814936765;

    REPt[0] = 4.56831690798003;   //Tension
    REPt[1] = 0.01374149543579;
	REPt[2] = 1.04617240539503;
	REPt[3] = 0.49323957503498;

	//Assign very small stiffness to slippage to avoid numerical problems
	Kslip = Keb / 10000;
	//u increment used for tangent evaluation
	u_inc = 0.0000001;

	// Variables subject to change
	// Initial Compressive Bearing Strain
	u_s0c = -Pslip / Keb;                     //compression slipping displacement
	f_s0c = -Pslip;

	u_b0c = -clearance + u_s0c;               //start of Richard Equation
	f_b0c = -clearance*Kslip + f_s0c;
	u_b0c_trial = u_b0c;
	f_b0c_trial = f_b0c;
	//Initial Tensile Bearing Strain
	u_s0t = Pslip / Keb;                      //tension slipping displacement
	f_s0t = Pslip;

	u_b0t = clearance + u_s0t;                //start of Richard Equation
	f_b0t = clearance*Kslip + f_s0t;
	u_b0t_trial = u_b0t;
	f_b0t_trial = f_b0t;
    // Internal Variables subjected to change
	u_bc = u_b0c;                           //start of bearing
	f_bc = f_b0c;
	u_bc_trial = u_bc;
	f_bc_trial = f_bc;

	u_bt = u_b0t;                           //start of bearing
	f_bt = f_b0t;
	u_bt_trial = u_bt;
	f_bt_trial = f_bt;

	u_sc = u_s0c;
	f_sc = f_s0c;
	u_st = u_s0t;
	f_st = f_s0t;

	//Initiate the breakpoints array
	ubkp[0]= u_b0c;
	ubkp[1] = u_b0c;
	ubkp[2] = u_s0c;
	ubkp[3] = u_s0t;
	ubkp[4] = u_b0t;
	ubkp[5] = u_b0t;
	
	fbkp[0] = f_b0c;
	fbkp[1] = f_b0c;
	fbkp[2] = f_s0c;
	fbkp[3] = f_s0t;
	fbkp[4] = f_b0t;
	fbkp[5] = f_b0t;

	//for u_current <= ubkp[1] and u_current >= ubkp[6] we use Richard Equation so there is no linear stiffness
	kbkp[0] = Kci;
	kbkp[1] = Kslip;
	kbkp[2] = Keb;
	kbkp[3] = Kslip;
	kbkp[4] = Kti;

	//status 1 tension bearing nonstick - 1 compression bearing nonstick
	//status 4 tension bearing stick - 4 compressiong bearing stick
	//status 0 before bearing
	status = 0;
	status_trial = status;
	//direction 1 loading - 1 reversal loading 0 for start
	dir = 0;
	dir_trial = dir;
	// CURRENT COMMITTED U  F  K
	u_current = 0.0;
	f_current = 0.0;
	k_current = E;
	// TRIAL U F K
	u_trial = 0.0;
	f_trial = 0.0;
	k_trial = E;

	u_transition = 0.0;
	f_transition = 0.0;
	u_transition_trial = u_transition;
	f_transition_trial = f_transition;

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();

}

SingleBoltAngleJoint::SingleBoltAngleJoint(void) :
	UniaxialMaterial(0, 0)
{
	//It's purpose is to create a blank  object, that will be filled in when the recvSelf() method is invoked on the object.
}

SingleBoltAngleJoint::~SingleBoltAngleJoint()
{
	// does nothing
}

UniaxialMaterial*
SingleBoltAngleJoint::getCopy(void)
{
	SingleBoltAngleJoint *theCopy = new SingleBoltAngleJoint(this->getTag(),
		Pslip,
		width_brace, thickness_brace, thickness_leg, Bolt_diameter, Hole_diameter, e1,
		Fy_brace, Fu_brace, Fy_leg,
		E, v);
	// FIXED PROPERTIES -----------------------------------------------------------------------------
	theCopy->clearance   = clearance; 
	theCopy->G           = G;
	theCopy->Leff        = Leff; 
	theCopy->An          = An; 
	theCopy->Aeff        = Aeff; 
	theCopy->Af          = Af; 
	theCopy->Anet        = Anet; 
	theCopy->Keb         = Keb; 
    // INITIAL COMPRESSIVE STIFFNESS
	theCopy->Kbr         = Kbr; 
	theCopy->Knet        = Knet;
	theCopy->Abolt       = Abolt; 
	theCopy->Lbolt       = Lbolt; 
	theCopy->Ibolt       = Ibolt;  
	theCopy->Kboltshear  = Kboltshear;
	theCopy->Kboltbending= Kboltbending; 
	theCopy->KbrLeg      = KbrLeg;
	theCopy->Kci         = Kci; 
	// INITIAL TENSILE STIFFNESS
	theCopy->Kbben       = Kbben;
	theCopy->Kbv         = Kbv; 
	theCopy->Kti         = Kti; 
	// COMPRESSIVE CAPACITY 
	theCopy->Rc          = Rc;
	// TENSILE CAPACITY
	theCopy->Rt          = Rt;
	// RICHARD EQUATION PARAMETER
	theCopy->REPc        = REPc;   
	theCopy->REPt        = REPt;   
    // OTHERS
	theCopy->Kslip       = Kslip; 
	theCopy->u_inc       = u_inc; 


    // VARIBLES SUBJECT TO CHANGE-------------------------------------------------------------------
    // INITIAL COMPRESSIVE BEARING DISPLACEMENT 
	theCopy->u_s0c       = u_s0c;  
	theCopy->f_s0c       = f_s0c;

	theCopy->u_b0c       = u_b0c;  
	theCopy->f_b0c       = f_b0c;
	theCopy->u_b0c_trial = u_b0c_trial;
	theCopy->f_b0c_trial = f_b0c_trial;
	// INITIAL TENSILE BEARING DISPLACEMENT 
	theCopy->u_s0t       = u_s0t;   
	theCopy->f_s0t       = f_s0t;

	theCopy->u_b0t       = u_b0t;  
	theCopy->f_b0t       = f_b0t;
	theCopy->u_b0t_trial = u_b0t_trial;
	theCopy->f_b0t_trial = f_b0t_trial;
	// bearing displacment and force
	theCopy->u_bc        = u_bc;
	theCopy->f_bc        = f_bc;
	theCopy->u_bc_trial = u_bc_trial;
	theCopy->f_bc_trial = f_bc_trial;

	theCopy->u_bt        = u_bt;
	theCopy->f_bt        = f_bt;
	theCopy->u_bt_trial = u_bt_trial;
	theCopy->f_bt_trial = f_bt_trial;
	// slipping displacement and force 
	theCopy->u_sc        = u_sc;
	theCopy->f_sc        = f_sc;
	theCopy->u_st        = u_st;
	theCopy->f_st        = f_st;

	// break points vector
	theCopy->ubkp        = ubkp;  
	theCopy->fbkp        = fbkp;  
	theCopy->kbkp        = kbkp;  

    // status 1 tension bearing nonstick -1 compression bearing nonstick 
    // status 4 tension bearing stick - 4 compressiong bearing stick
    // status 0 before bearing
	theCopy->status      = status;
	theCopy->status_trial = status_trial;

	//direction 1 loading - 1 reversal loading 0 for start
	theCopy->dir         = dir;
	theCopy->dir_trial = dir_trial;

	// CURRENT COMMITTED U  F  K
	theCopy->u_current   = u_current;
	theCopy->f_current   = f_current;
	theCopy->k_current   = k_current; 
	// TRIAL U F K
	theCopy->u_trial     = u_trial; // can be defined as trialStrain from setTrailStrain
	theCopy->f_trial     = f_trial;
	theCopy->k_trial     = k_trial;
	// TRANSITION U F FOR BEARING STICK PHASE
	theCopy->u_transition= u_transition;
	theCopy->f_transition= f_transition;
	theCopy->u_transition_trial = u_transition_trial;
	theCopy->f_transition_trial = f_transition_trial;


	return theCopy;
}


//functions
double SingleBoltAngleJoint::Line(double x, double k, double u, double f)
{
	//calculates the y on line (k, u, f) with x
	return  k*x + (f - k*u);

}

double SingleBoltAngleJoint::Intersectu(double k1, double u1, double f1, double k2, double u2, double f2)
{
	//calculates the deformation of intersection point between two lines 
	return ((f2 - k2*u2) - (f1 - k1*u1)) / (k1 - k2);

}

double SingleBoltAngleJoint::Intersectf(double k1, double u1, double f1, double k2, double u2, double f2)
{
	//calculates the force of intersection point between two lines 
	return k1*((f2 - k2*u2) - (f1 - k1*u1)) / (k1 - k2) + (f1 - k1*u1);

}

double SingleBoltAngleJoint::Bisectionu(double u, double f, double Kb, double_vec REP, double u_b0, double R, double Ki, int status, double f_b0, double Kslip)
{
	//Return the displacement of intersection of lin (u, f, Ki) and RichardEquation 
	double low;
	double high;
	double guess;
	double maxIterations;
	double numIterations;
	double F1;
	double F2;
	double error;
	double tol = 0.00001;
	if (u > 0) {
		low = u;   //tension
		high = u + 1;
	}
	else {
		low = u;   //compression
		high = u - 1;
	}
	guess = (low + high) / 2;
	maxIterations = 100000;
	numIterations = 0;
	while (numIterations < maxIterations) {
		F1 = RichardEquation(REP, f_b0, Kslip, R, Ki, guess, u_b0);

		if (F1 < 0)
		{   //compression
			F2 = f + (guess - u)*Kb;
		}
		if (F1 > 0)
		{// tension
			F2 = f + (guess - u)*Kb;
		}
		error = abs(F2) - abs(F1);
		if (abs(error) < tol)
		{
			break;
		}

		if (error > 0)
		{
			high = guess;
		}
		else
		{
			low = guess;
		}
		guess = (low + high) / 2;
		numIterations = numIterations + 1;
	}
	
	return guess;

}


double SingleBoltAngleJoint::Bisectionf(double u, double f, double Kb, double bisectionu)
{ //Return the force of lin (u, f, Ki) and RichardEquation 

	return f + (bisectionu - u)*Kb;

}

double SingleBoltAngleJoint::RichardEquation(double_vec REP, double f_b0, double Kslip, double R, double Ki, double u, double u_b0)
{ //Compute the stress based on given Richard Equation 
	double def_norm;
	double a;
	double b;
	double F;
	def_norm = abs((u - u_b0)*Ki / (R - abs(f_b0)));  //positive
	a = pow (1 + def_norm*REP[0] / REP[2], REP[3]);
	b = pow(a, 1 / REP[3]);
	F= (def_norm*REP[0]/b+ def_norm*REP[1])*(R - abs(f_b0)) + abs(f_b0);
	if (f_b0 > 0)
	{
		F = F + Kslip*(u - u_b0);   //positive for tension
	}
	else
	{
		F = -F + Kslip*(u - u_b0);  //negative for compression
	}

	return F;
}


void SingleBoltAngleJoint::CalcBreakpoints(int dir, int status, double Pslip, double u_bc, double f_bc, double u_bt, double f_bt, double u_b0c, double f_b0c,
	double u_b0t, double f_b0t, double u_s0c, double f_s0c, double u_s0t, double f_s0t, double u_current, double f_current,
	double Keb, double Kti, double Kci, double_vec REPc, double_vec REPt, double Rc, double Rt, double Kslip, double u_transition, double f_transition,
	double_vec& ubkp, double_vec& fbkp, double_vec& kbkp)


{  // compute control points for force calculation 
	double tol = 0.0000001;
	// Case 1 u_current before bearing
	if ((f_current <= Line(u_current, Keb, u_bc, f_bc)) && (f_current >= Line(u_current, Keb, u_bt, f_bt)))
	{
		if (u_bc == u_b0c)
		{ //there is no  enlongation in the compression yet
			ubkp[0] = u_bc;
			fbkp[0] = f_bc;
		}
		else
		{
			ubkp[0] = Bisectionu(u_bc, f_bc, Kci, REPc, u_b0c, Rc, Kci, status, f_b0c, Kslip);  //u_bc intersection with RE at Compression
			fbkp[0] = Bisectionf(u_bc, f_bc, Kci, ubkp[0]);
		}
		ubkp[1] = u_bc;
		fbkp[1] = f_bc;
		ubkp[2] = Intersectu(Kslip, u_s0c, f_s0c, Keb, u_current, f_current); //intersect between compression slipping and u_current
		fbkp[2] = Intersectf(Kslip, u_s0c, f_s0c, Keb, u_current, f_current); //intersect between compression slipping and u_current
		ubkp[3] = Intersectu(Kslip, u_s0t, f_s0t, Keb, u_current, f_current); //intersect between compression slipping and u_current
		fbkp[3] = Intersectf(Kslip, u_s0t, f_s0t, Keb, u_current, f_current); //intersect between compression slipping and u_current
		ubkp[4] = u_bt;
		fbkp[4] = f_bt;
		if (u_bt == u_b0t)
		{   //there is no  enlongation in the tension yet
			ubkp[5] = u_bt;
			fbkp[5] = f_bt;
		}
		else
		{
			ubkp[5] = Bisectionu(u_bt, f_bt, Kti, REPt, u_b0t, Rt, Kti, status, f_b0t, Kslip);  //u_bt intersection with RE at Tension
			fbkp[5] = Bisectionf(u_bt, f_bt, Kti, ubkp[5]);
		}
		kbkp[0] = Kci;
		kbkp[1] = Kslip;
		kbkp[2] = Keb;
		kbkp[3] = Kslip;
		kbkp[4] = Kti;
	}

	//Case 2 u_current after tension bearing 
	else if (status > 0)
	{

		if (u_bc == u_b0c) //there is no  enlongation in the compression yet
		{
			ubkp[0] = u_bc;
			fbkp[0] = f_bc;
		}
		else
		{
			ubkp[0] = Bisectionu(u_bc, f_bc, Kci, REPc, u_b0c, Rc, Kci, status, f_b0c, Kslip);  //u_bc intersection with RE at Compression
			fbkp[0] = Bisectionf(u_bc, f_bc, Kci, ubkp[0]);
		}

		ubkp[1] = u_bc;
		fbkp[1] = f_bc;
		if (status == 4)  // stick
		{
			if (f_transition >= f_current)
			{
				ubkp[4] = u_transition;                               // equals to transition point
				fbkp[4] = f_transition;
				ubkp[3] = u_transition - 2 * Pslip / Keb;             //transition point 2Pslip downward with Keb
				fbkp[3] = f_transition - 2 * Pslip;
			}
			else
			{
				ubkp[3] = u_transition;                               // equals to transition point
				fbkp[3] = f_transition;
				ubkp[4] = u_transition + 2 * Pslip / Keb;             // transition point 2Pslip upward with Keb
				fbkp[4] = f_transition + 2 * Pslip;
			}
		}

		else if (status != 4)
		{// non stick
			if (dir == 1)
			{// loading
				ubkp[3] = u_current - 2 * Pslip / Keb;             // transition point 2Pslip downward with Keb
				fbkp[3] = f_current - 2 * Pslip;
				ubkp[4] = u_current;                               // equals to current point
				fbkp[4] = f_current;
			}
			else if (dir == -1)
			{
				ubkp[3] = u_current;                               // equals to transition point
				fbkp[3] = f_current;
				ubkp[4] = u_current + 2 * Pslip / Keb;             // transition point 2Pslip upward with Keb
				fbkp[4] = f_current + 2 * Pslip;
			}
		}
		if (abs(fbkp[4] - RichardEquation(REPt, f_b0t, Kslip, Rt, Kti, ubkp[4], u_b0t)) < tol)  //break point 5 equals and larger than break point 6
		{
			ubkp[5] = ubkp[4];
			fbkp[5] = fbkp[4];
		}
		else
		{
			ubkp[5] = Bisectionu(ubkp[4], fbkp[4], Kti, REPt, u_b0t, Rt, Kti, status, f_b0t, Kslip);  //u_bt intersection with RE at Tension
			fbkp[5] = Bisectionf(ubkp[4], fbkp[4], Kti, ubkp[5]);
		}
		ubkp[2] = Intersectu(Kslip, u_s0c, f_s0c, Kti, ubkp[3], fbkp[3]); //intersect between compression slipping and point4
		fbkp[2] = Intersectf(Kslip, u_s0c, f_s0c, Kti, ubkp[3], fbkp[3]); //intersect between compression slipping and point4
		kbkp[0] = Kci;
		kbkp[1] = Kslip;
		kbkp[2] = Kti;
		kbkp[3] = Keb;
		kbkp[4] = Kti;
	}

	//Case 3 u_current after compression bearing 
	else if (status < 0)
	{
		if (u_bt == u_b0t)
		{//there is no  enlongation in the tension yet
			ubkp[5] = u_bt;
			fbkp[5] = f_bt;
		}
		else
		{
			ubkp[5] = Bisectionu(u_bt, f_bt, Kti, REPt, u_b0t, Rt, Kti, status, f_b0t, Kslip);  //u_bt intersection with RE at Tension
			fbkp[5] = Bisectionf(u_bt, f_bt, Kti, ubkp[5]);
		}
		ubkp[4] = u_bt;
		fbkp[4] = f_bt;

		if (status == -4)
		{
			if (f_transition >= f_current)
			{
				ubkp[2] = u_transition;            // equals to transition point
				fbkp[2] = f_transition;
				ubkp[1] = u_transition - 2 * Pslip / Keb;             // transition point 2Pslip downward with Keb
				fbkp[1] = f_transition - 2 * Pslip;
			}
			else
			{
				ubkp[1] = u_transition;            // equals to transition point
				fbkp[1] = f_transition;
				ubkp[2] = u_transition + 2 * Pslip / Keb;             // transition point 2Pslip upward with Keb
				fbkp[2] = f_transition + 2 * Pslip;
			}
		}

		else if (status != -4)
		{
			if (dir == 1)
			{// loading
				ubkp[1] = u_current;            // equals to current point
				fbkp[1] = f_current;
				ubkp[2] = u_current + 2 * Pslip / Keb;             // current point 2Pslip downward with Keb
				fbkp[2] = f_current + 2 * Pslip;
			}
			else if (dir == -1)
			{ // unloading
				ubkp[2] = u_current;            // equals to current point
				fbkp[2] = f_current;
				ubkp[1] = u_current - 2 * Pslip / Keb;             // transition point 2Pslip downward with Keb
				fbkp[1] = f_current - 2 * Pslip;
			}
		}

		if (abs(fbkp[1] - RichardEquation(REPc, f_b0c, Kslip, Rc, Kci, ubkp[1], u_b0c)) < tol)
		{//break point 2 equals and smaller than break point 1
			ubkp[0] = ubkp[1];
			fbkp[0] = fbkp[1];
		}
		else
		{
			ubkp[0] = Bisectionu(ubkp[1], fbkp[1], Kci, REPc, u_b0c, Rc, Kci, status, f_b0c, Kslip);  //u_bc intersection with RE at Compression
			fbkp[0] = Bisectionf(ubkp[1], fbkp[1], Kci, ubkp[0]);
		}
		ubkp[3] = Intersectu(Kslip, u_s0t, f_s0t, Kci, ubkp[2], fbkp[2]); //intersect between tension slipping and point
		fbkp[3] = Intersectf(Kslip, u_s0t, f_s0t, Kci, ubkp[2], fbkp[2]); //intersect between tension slipping and point
		kbkp[0] = Kci;
		kbkp[1] = Keb;
		kbkp[2] = Kci;
		kbkp[3] = Kslip;
		kbkp[4] = Kti;
	}
}

void SingleBoltAngleJoint::CalcForceStatusDir(double_vec ubkp, double_vec fbkp, double_vec kbkp, double u_trial, int status, int dir, double u_current, double f_current, double u_st, double f_st,
	double u_sc, double f_sc, double Kslip, double Kti, double Kci, double u_transition, double f_transition,
	double u_bc, double f_bc, double u_bt, double f_bt, double_vec REPt, double_vec REPc, double Rt, double Rc, double u_b0t, double f_b0t, double u_b0c, double f_b0c, double Keb, double Pslip,
	double& f_trial, int& status_trial, int& dir_trial, double& u_transition_trial, double& f_transition_trial, double& u_bc_trial, double& f_bc_trial, double& u_bt_trial,
	double& f_bt_trial, double& u_b0t_trial, double& f_b0t_trial, double& u_b0c_trial, double& f_b0c_trial, double& k_trial)
{
	u_transition_trial = u_transition;
	f_transition_trial = f_transition;
	u_bc_trial = u_bc;
	f_bc_trial = f_bc;
	u_bt_trial = u_bt;
	f_bt_trial = f_bt;
	u_b0t_trial = u_b0t;
	f_b0t_trial = f_b0t;
	u_b0c_trial = u_b0c;
	f_b0c_trial = f_b0c;


	if (u_trial < ubkp[0])
	{
		f_trial = RichardEquation(REPc, f_b0c, Kslip, Rc, Kci, u_trial, u_b0c);
	}
	else if ((u_trial >= ubkp[0]) && (u_trial < ubkp[1]))
	{
		f_trial = fbkp[1] + kbkp[0]*(u_trial - ubkp[1]);
	}
	else if ((u_trial >= ubkp[1]) && (u_trial < ubkp[2]))
	{
		f_trial = fbkp[2] + kbkp[1]*(u_trial - ubkp[2]);
	}
	else if ((u_trial >= ubkp[2]) && (u_trial < ubkp[3]))
	{
		f_trial = fbkp[3] + kbkp[2]*(u_trial - ubkp[3]);
	}
	else if ((u_trial >= ubkp[3]) && (u_trial < ubkp[4]))
	{
		f_trial = fbkp[4] + kbkp[3]*(u_trial - ubkp[4]);
	}
	else if ((u_trial >= ubkp[4]) && (u_trial < ubkp[5]))
	{
		f_trial = fbkp[5] + kbkp[4]*(u_trial - ubkp[5]);
	}
	else
	{
		f_trial = RichardEquation(REPt, f_b0t, Kslip, Rt, Kti, u_trial, u_b0t);
	}

	switch (status) {

	case 0: // current point not in bearing stage
		if (u_trial > ubkp[4])
		{
			status_trial = 1;
		}
		else if (u_trial < ubkp[1])
		{
			status_trial = -1;
		}
		else
		{
			status_trial = status; //status not change
		}

	case 1: // current point in the bearing nonstick tension
		if (u_trial >= ubkp[4])
		{
			status_trial = 1;
		}
		else if ((u_trial < ubkp[4]) && (u_trial >= ubkp[3]))
		{
			status_trial = 4;
			if (dir == 1)
			{
				u_transition_trial = ubkp[4];
				f_transition_trial = fbkp[4];
			}
			else if (dir == -1)
			{
				u_transition_trial = ubkp[3];
				f_transition_trial = fbkp[3];
			}
		}
		else if ((u_trial < ubkp[3]) && (u_trial >= ubkp[2]))
		{
			status_trial = 1;
		}
		else if ((u_trial < ubkp[2]) && (u_trial >= ubkp[1]))
		{
			status_trial = 0;
			// update the new tension bearing point
			u_bt_trial = Intersectu(Keb, ubkp[2], fbkp[2], Kslip, u_st, f_st);
			f_bt_trial = Intersectu(Keb, ubkp[2], fbkp[2], Kslip, u_st, f_st);
			//comporession bearing point upadate !!!!!!!!!!!!!!!!!!!!!!!
			u_bc_trial = (u_bt_trial - u_bt)*0.2 + u_bc; //20 % is not going to elongation
			f_bc_trial = f_sc + Kslip*(u_bc_trial - u_sc);
			u_b0c_trial = (u_bt_trial - u_bt)*0.2 + u_b0c;
			f_b0c_trial = f_sc + Kslip*(u_b0c_trial - u_sc);
		}
		else if (u_trial < ubkp[1])
		{
			status_trial = -1;
		}

	case 4: // current point in the bearing stick tension
		if (u_trial >= ubkp[4])
		{
			status_trial = 1;
		}
		else if ((u_trial < ubkp[4]) && (u_trial >= ubkp[3]))
		{
		    status_trial = 4;
		}
		else if ((u_trial < ubkp[3]) && (u_trial >= ubkp[2]))
		{
			status_trial = 1;
		}
		else if ((u_trial < ubkp[2]) && (u_trial >= ubkp[1]))
		{
			status_trial = 0;
		}
		else if (u_trial < ubkp[1])
		{
			status_trial = -1;
		}
	case -1: //current point in the bearing nonstick compression
		if (u_trial >= ubkp[4])
		{
			status_trial = 1;
		}
		else if ((u_trial < ubkp[4]) && (u_trial >= ubkp[3]))
		{
			status_trial = 0;
			//update the new compression bearing point
			u_bc_trial = Intersectu(Keb, ubkp[2], fbkp[2], Kslip, u_sc, f_sc);
			f_bc_trial = Intersectf(Keb, ubkp[2], fbkp[2], Kslip, u_sc, f_sc);
			//tension bearing point upadate !!!!!!!!!!!!!!!!!!!!!!!
			u_bt_trial = (u_bc_trial - u_bc)*0.4 + u_bt; // 40 % is not going to elongation
			f_bt_trial = f_st + Kslip*(u_bt_trial - u_st);
			u_b0t_trial = (u_bc_trial - u_bc)*0.4 + u_b0t;
			f_b0t_trial = f_st + Kslip*(u_b0t_trial - u_st);
		}
		else if ((u_trial < ubkp[3]) && (u_trial >= ubkp[2]))
		{
			status_trial = -1;
		}
		else if ((u_trial < ubkp[2]) && (u_trial >= ubkp[1]))
		{
			status_trial = -4;
			if (dir == 1)
			{
				u_transition_trial = ubkp[1];
				f_transition_trial = fbkp[1];
			}
			else if (dir == -1)
			{
				u_transition_trial = ubkp[2];
				f_transition_trial = fbkp[2];
			}
		}
		else if (u_trial < ubkp[1])
		{
			status_trial = -1;
		}

	case -4: // current point in the bearing stick tension
		if (u_trial >= ubkp[4])
		{
			status_trial = 1;
		}
		else if ((u_trial < ubkp[4]) && (u_trial >= ubkp[3]))
		{
			status_trial = 0;
		}
		else if ((u_trial < ubkp[3]) && (u_trial >= ubkp[2]))
		{
			status_trial = -1;
		}
		else if ((u_trial < ubkp[2]) && (u_trial >= ubkp[1]))
		{
			status_trial = -4;
		}
		else if (u_trial < ubkp[1])
		{
			status_trial = -1;
		}



	}

	if ((u_trial >= u_bt_trial && f_trial > 0) || (u_trial >= u_bt_trial - 2 * Pslip / Keb && f_trial < 0)) // in tension bearing
	{
		if (f_trial - f_current >= 0)
		{
			dir_trial = 1;
		}
		else
		{
			dir_trial = -1;
		}
	}
	else if ((u_trial <= u_bc_trial && f_trial < 0) || (u_trial <= u_bc_trial + 2 * Pslip / Keb && f_trial > 0))
	{//in compression bearing
		if (f_trial - f_current >= 0)
		{
			dir_trial = -1;
		}
		else
		{
			dir_trial = 1;
		}
	}
	else
	{
		dir_trial = 0;
	}
	// compute the stress increase over u_inc for stiffness assessment 
	double f_stiff;
	double u_stiff = u_trial - u_inc;
	if (u_stiff < ubkp[0])
	{
		f_stiff = RichardEquation(REPc, f_b0c, Kslip, Rc, Kci, u_stiff, u_b0c);
	}
	else if ((u_stiff >= ubkp[0]) && (u_stiff < ubkp[1]))
	{
		f_stiff = fbkp[1] + kbkp[0] * (u_stiff - ubkp[1]);
	}
	else if ((u_stiff >= ubkp[1]) && (u_stiff < ubkp[2]))
	{
		f_stiff = fbkp[2] + kbkp[1] * (u_stiff - ubkp[2]);
	}
	else if ((u_stiff >= ubkp[2]) && (u_stiff < ubkp[3]))
	{
		f_stiff = fbkp[3] + kbkp[2] * (u_stiff - ubkp[3]);
	}
	else if ((u_stiff >= ubkp[3]) && (u_stiff < ubkp[4]))
	{
		f_stiff = fbkp[4] + kbkp[3] * (u_stiff - ubkp[4]);
	}
	else if ((u_stiff >= ubkp[4]) && (u_stiff < ubkp[5]))
	{
		f_stiff = fbkp[5] + kbkp[4] * (u_stiff - ubkp[5]);
	}
	else
	{
		f_stiff = RichardEquation(REPt, f_b0t, Kslip, Rt, Kti, u_stiff, u_b0t);
	}
	k_trial = (f_trial - f_stiff) / u_inc;

}

int SingleBoltAngleJoint::setTrialStrain(double strain, double strainRate )
{  
   //if delta u is very small, there is nothing to calculate
	if (fabs(u_trial - strain) < DBL_EPSILON)
	{
		return 0;
	}
	// assign the new displacement 
	u_trial = strain;

	//calculate the break points and update
	CalcBreakpoints( dir,  status,  Pslip,  u_bc,  f_bc,  u_bt,  f_bt,  u_b0c,  f_b0c,
		 u_b0t,  f_b0t,  u_s0c,  f_s0c,  u_s0t,  f_s0t,  u_current,  f_current,
		 Keb,  Kti,  Kci,  REPc,  REPt,  Rc,  Rt,  Kslip,  u_transition,  f_transition,
		 ubkp, fbkp,  kbkp);
	//caluate the trial variables for force status dir and etc. update
	CalcForceStatusDir( ubkp,  fbkp, kbkp,  u_trial,  status,  dir,  u_current,  f_current,  u_st,  f_st,
		 u_sc,  f_sc,  Kslip,  Kti,  Kci,  u_transition,  f_transition,
		 u_bc,  f_bc,  u_bt,  f_bt,  REPt,  REPc,  Rt,  Rc,  u_b0t,  f_b0t,  u_b0c,  f_b0c,  Keb,  Pslip,
		 f_trial,  status_trial,  dir_trial,  u_transition_trial,  f_transition_trial,  u_bc_trial,  f_bc_trial,  u_bt_trial,
		 f_bt_trial,  u_b0t_trial,  f_b0t_trial,  u_b0c_trial,  f_b0c_trial, k_trial);


	return 0;
}

double SingleBoltAngleJoint::getStrain(void)
{
	return u_trial;
}
double SingleBoltAngleJoint::getStress(void)
{
	return f_trial;
}
double SingleBoltAngleJoint::getTangent(void)
{
	return k_trial;
}
int SingleBoltAngleJoint::commitState(void)
{
	u_current = u_trial;
	f_current = f_trial;
	k_current = k_trial;
	// store the variables subject to change 
	status = status_trial;
	dir = dir_trial;
	u_transition = u_transition_trial;
	f_transition = f_transition_trial;
	u_bc = u_bc_trial;
	f_bc = f_bc_trial;
	u_bt = u_bt_trial;
	f_bt = f_bt_trial;
	u_b0t = u_b0t_trial;
	f_b0t = f_b0t_trial;
	u_b0c = u_b0c_trial;
	f_b0c = f_b0c_trial;

	return 0;
}

double SingleBoltAngleJoint::getInitialTangent(void)
{
	return E;
}


int SingleBoltAngleJoint::revertToLastCommit(void)
{   // Load the stored values
	u_trial = u_current ;
	f_trial = f_current ;
	k_trial = k_current ;

	status_trial = status;
	dir_trial = dir;
	u_transition_trial = u_transition;
	f_transition_trial = f_transition;
	u_bc_trial = u_bc;
	f_bc_trial = f_bc;
	u_bt_trial = u_bt;
	f_bt_trial = f_bt;
	u_b0t_trial = u_b0t;
	f_b0t_trial = f_b0t;
	u_b0c_trial = u_b0c;
	f_b0c_trial = f_b0c;

	return 0;
}

int SingleBoltAngleJoint::revertToStart(void)
{   // initialize the stored values
	u_trial = u_current = 0.0;
	f_trial = f_current = 0.0;
	k_trial = k_current = E;

	status_trial = status = 0;
	dir_trial = dir = 0;
	u_transition_trial = u_transition = 0.0;
	f_transition_trial = f_transition = 0.0;
	u_bc_trial = u_bc = -clearance + u_s0c;
	f_bc_trial = f_bc = -clearance*Kslip + f_s0c;
	u_bt_trial = u_bt = clearance + u_s0t;
	f_bt_trial = f_bt = clearance*Kslip + f_s0t;
	u_b0t_trial = u_b0t = clearance + u_s0t;
	f_b0t_trial = f_b0t = clearance*Kslip + f_s0t;
	u_b0c_trial = u_b0c = -clearance + u_s0c;
	f_b0c_trial = f_b0c = -clearance*Kslip + f_s0c;

	return 0;
}

int SingleBoltAngleJoint::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int SingleBoltAngleJoint::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	return -1;
}
void SingleBoltAngleJoint::Print(OPS_Stream &s, int flag )
{
	s << "SingleBoltAngleJoint tag: " << this->getTag() << endln;
	s << "  deformation: " << u_trial << "  force: " << f_trial << " tangent: " <<k_trial << endln;

}













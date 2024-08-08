/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

/*
rotate / pitch body depending on torque moment on z axis
y vertical translation
2dof
solid body rotation, no deformation
uses newmark scheme and internal forces and moments
no reading of forces, moments from external file
added saved_time_
*/

#include "rotatingMotion_3dof_pitch_xy_translate_ext_thrust.H"
#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	namespace solidBodyMotionFunctions
	{
		defineTypeNameAndDebug(rotatingMotion_3dof_pitch_xy_translate_ext_thrust, 0);
		addToRunTimeSelectionTable
		(
			solidBodyMotionFunction,
			rotatingMotion_3dof_pitch_xy_translate_ext_thrust,
			dictionary
		);
	}
	
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingMotion_3dof_pitch_xy_translate_ext_thrust::rotatingMotion_3dof_pitch_xy_translate_ext_thrust
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
	dict_(SBMFCoeffs),
    //motion_(SBMFCoeffs,SBMFCoeffs),
	motion_(SBMFCoeffs,SBMFCoeffs,time_.time()),
	firstRun_(true),
    //startTimeIndex_(-1),
	//waitToRunIndex_(SBMFCoeffs_.get<label>("waitToRunIndex")),
	//waitToRunIndexRamp_(SBMFCoeffs_.get<label>("waitToRunIndexRamp")),
	//fileWriteInterval_(SBMFCoeffs_.get<scalar>("fileWriteInterval")),
	fileWriteInterval_(runTime.controlDict().get<scalar>("writeInterval")),
	angle_old(0),
	angular_vel_z_old(0),
	angular_acc_z_old(0),
	//angle_ramp(0),
	//angular_vel_z_ramp(0),
	//angular_acc_z_ramp(0),
	patches_(wordRes(SBMFCoeffs.lookup("patches"))),
	rhoInf_(SBMFCoeffs_.getOrDefault<scalar>("rhoInf", 1.0)),
	rhoName_(SBMFCoeffs.lookupOrDefault<word>("rho", "rho")),
    origin_(SBMFCoeffs_.get<vector>("origin")),
    axis_(SBMFCoeffs_.get<vector>("axis")),
	mass_(SBMFCoeffs_.get<scalar>("mass")),
	centreOfMass_(SBMFCoeffs_.get<vector>("centreOfMass")),
    momentOfInertia_(SBMFCoeffs_.get<vector>("momentOfInertia")),
    accelerationRelaxation_(SBMFCoeffs_.getOrDefault<scalar>("accelerationRelaxation", 1.0)),
	accelerationDamping_(SBMFCoeffs_.getOrDefault<scalar>("accelerationDamping", 1.0)),
	g_(SBMFCoeffs_.get<vector>("g")),
	gamma_new_(SBMFCoeffs_.getOrDefault<scalar>("gamma_new", 0.5)),
	beta_new_
	(
		max
		(
		 0.25*(gamma_new_ + 0.5)*(gamma_new_ + 0.5),
		 SBMFCoeffs_.getOrDefault<scalar>("beta_new", 0.25)
		)
	),
	initial_velocity_(SBMFCoeffs_.get<vector>("initial_velocity")),
	x_displ_old(0),
	vel_x_old(0),
	acc_x_old(0),
	y_displ_old(0),
	vel_y_old(0),
	acc_y_old(0),
	
	saved_time_(-1),
	saved_no(1),
	angleEveryStep_(SBMFCoeffs_.getOrDefault<bool>("angleEveryStep", 0)),
	max_propeller_thrust_(SBMFCoeffs_.getOrDefault<scalar>("max_propeller_thrust", 1.0)),
	activation_pitch_angle_(SBMFCoeffs_.getOrDefault<scalar>("activation_pitch_angle", 1.0)),
	activation_time_(SBMFCoeffs_.getOrDefault<scalar>("activation_time", 1.0)),
	propeller_coeff_(SBMFCoeffs_.getOrDefault<scalar>("propeller_coeff", 1.0))
    //omega_(Function1<scalar>::New("omega", SBMFCoeffs_, &runTime))
{
    read(SBMFCoeffs);

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotatingMotion_3dof_pitch_xy_translate_ext_thrust::transformation() const
{
	
	string line;
	scalar tmp_value;
	//scalar tmp_value1;
	//scalar tmp_value2;
	//scalar tmp_value3;
	scalar angular_acc_z;
	//scalar angular_acc_z_old;
	scalar angular_vel_z;
	//scalar angular_vel_z_old;
	scalar angle;
	//scalar angle_old;
	scalar delta_time;
	scalar acc_x;
	scalar vel_x;
	scalar x_displ;
	scalar delta_x;
	scalar acc_y;
	scalar vel_y;
	scalar y_displ;
	scalar delta_y;
    scalar t = time_.value();
	//scalar write_interval = time_.writeInterval();
	vector displacement;
		
	if (firstRun_ == true) {
	
		startTimeIndex_ = time_.timeIndex(); //+ waitToRunIndex_;
		
		firstRun_ = false;
		
		tmp_value = 0;
		//tmp_value1 = 0;
		//tmp_value2 = 0;
		//tmp_value3 = 0;
		angular_acc_z = 0;
		//angular_acc_z_old = 0;
		angular_vel_z = 0;
		//angular_vel_z_old = 0;
		angle = 0;
		//angle_old = 0;
		delta_time = 0;
		saved_time_ = -1;	//just any negative value
		saved_no = t/fileWriteInterval_ + 1;
		
		acc_x = 0;
		vel_x = 0;
		x_displ = 0;
		
		acc_y = 0;
		vel_y = 0;
		y_displ = 0;
		displacement = vector(0,0,0);
		
		//origin_ = vector(0,0,0);
		//centreOfMass_ = vector(0,0,0);
		
		//Info << "startTimeIndex = " << startTimeIndex_ << " fileWriteInterval = " << fileWriteInterval_ << " waitToRunIndex_ = " << waitToRunIndex_ << endl;
		
		Info << "patches = " << patches_ << " rhoInf = " << rhoInf_ << " CofR = " << motion_.centreOfRotation() << " rho = " << rhoName_ << endl;
				
		Info << "mass = " << mass_ << " momentOfInertia = " << momentOfInertia_[2] << endl;
		
		Info << "gamma_new = " << gamma_new_ << " beta_new = " << beta_new_ << " g = " << g_ << endl;
		
		Info << "accelerationRelaxation = " << accelerationRelaxation_ << " accelerationDamping = " << accelerationDamping_ << endl;

		Info << "activation_pitch_angle = " << activation_pitch_angle_ << " activation_time = " << activation_time_ << " propeller_coef = " << propeller_coeff_ << endl;
		
	}
	
	
	//Info << "startTimeIndex = " << startTimeIndex_ << " time_.timeIndex() = " << time_.timeIndex() << " waitToRunIndex_ = " << waitToRunIndex_ << endl;
	//Info << "fileWriteInterval = " << fileWriteInterval_ << endl;
	
	if (t - saved_time_ < 1e-8 && time_.timeIndex() != 0) {
		
		angle = angle_old;
		angular_vel_z = angular_vel_z_old;
		angular_acc_z = angular_acc_z_old;
		
		x_displ = x_displ_old;
		vel_x = vel_x_old;
		acc_x = acc_x_old;
		y_displ = y_displ_old;
		vel_y = vel_y_old;
		acc_y = acc_y_old;
		
	}
		
	else {
		
		Info << "rotatingMotion_3dof_pitch_xy_translate_ext_thrust time_.timeIndex() = " << time_.timeIndex() << endl;
	
		//normal run, except when time_.timeIndex() == startTimeIndex_
		//if (time_.timeIndex() >= startTimeIndex_ ) {
		
		
		// add force
		dictionary forcesDict;
		forcesDict.add("type", functionObjects::forces::typeName);
		forcesDict.add("patches", patches_);

		forcesDict.add("rhoInf", rhoInf_);
		forcesDict.add("rho", rhoName_);
		forcesDict.add("CofR", motion_.centreOfRotation());


		functionObjects::forces f("forces", time_, forcesDict);
		
		f.calcForcesMoments();
		Info << "force = " << f.forceEff() << " moment = " << f.momentEff() << endl;
		//Info << "patches = " << patches_ << " rhoInf = " << rhoInf_ << " CofR = " << motion_.centreOfRotation() << " rho = " << rhoName_ << endl;
			
		delta_time = time_.deltaTValue();
		
		//Info << "mass = " << mass_ << " momentOfInertia = " << momentOfInertia_[2] << " f.momentEff()[2] = " << f.momentEff()[2] << endl;
		
		scalar propeller_thrust_x = 0;
		scalar propeller_thrust_y = 0;
		scalar propeller_thrust = 0;
		
		if (time_.timeIndex() != 1) {
			
			if (time_.timeIndex() != startTimeIndex_) {
			
				Info << "6dof applied with newmark" << endl;
				
				//read in saved time and incoming velocity
				/*ifstream MyReadFile("current_time_angle_vel_acc.dat");

				while(MyReadFile.good() && (getline(MyReadFile, line)))
				{
					istringstream iss(line);

					iss >> std::fixed >> std::setprecision(8) >> tmp_value >> angle_old >> angular_vel_z_old >> angular_acc_z_old;

				}	

				MyReadFile.close();*/
				
				//Info << "angle_old = " << angle_old << " angular_vel_z_old = " << angular_vel_z_old << " angular_acc_z_old = " << angular_acc_z_old << endl;
				
				//Info << "y_displ_old = " << y_displ_old << " vel_y_old = " << vel_y_old << " acc_y_old = " << acc_y_old << endl;

			}
			
			else {
				
				//just resumed
					
				Info << "resumed using prev angle, displ, vel, acc. together with 6dof" << endl;
				
				//read in saved time and incoming velocity
				ifstream MyReadFile("interval_angle_displ_vel_acc_CG.dat");

				while(MyReadFile.good() && (getline(MyReadFile, line))) {
					
					istringstream iss(line);

					iss >> std::fixed >> std::setprecision(8) >> tmp_value >> angle_old >> angular_vel_z_old >> angular_acc_z_old >> x_displ_old >> vel_x_old >> acc_x_old >> y_displ_old >> vel_y_old >> acc_y_old >> motion_.centreOfRotation().component(vector::X) >> motion_.centreOfRotation().component(vector::Y) >> motion_.centreOfRotation().component(vector::Z);

				}
				
				centreOfMass_ = motion_.centreOfRotation();
				
			}
			
			Info << "angle_old = " << angle_old << " angular_vel_z_old = " << angular_vel_z_old << " angular_acc_z_old = " << angular_acc_z_old << endl;
			
			Info << "x_displ_old = " << x_displ_old << " vel_x_old = " << vel_x_old << " acc_x_old = " << acc_x_old << endl;
			
			Info << "y_displ_old = " << y_displ_old << " vel_y_old = " << vel_y_old << " acc_y_old = " << acc_y_old << endl;
			
			Info << "motion_.centreOfRotation() = " << motion_.centreOfRotation() << endl;

			angular_acc_z = accelerationRelaxation_*(f.momentEff()[2]/momentOfInertia_[2]) + (1. - accelerationRelaxation_)*angular_acc_z_old;
			
			angular_vel_z = angular_vel_z_old + accelerationDamping_*delta_time*(gamma_new_*angular_acc_z + (1. - gamma_new_)*angular_acc_z_old);
			
			angle = angle_old + angular_vel_z_old*delta_time + accelerationDamping_*delta_time*delta_time*(beta_new_*angular_acc_z + (0.5 - beta_new_)*angular_acc_z_old);
			
			acc_x = accelerationRelaxation_*(f.forceEff()[0]/mass_ + g_[0]) + (1. - accelerationRelaxation_)*acc_x_old;
			
			vel_x = vel_x_old + accelerationDamping_*delta_time*(gamma_new_*acc_x + (1. - gamma_new_)*acc_x_old);
			
			delta_x = vel_x*delta_time + accelerationDamping_*delta_time*delta_time*(beta_new_*acc_x + (0.5 - beta_new_)*acc_x_old);
			
			
			
			//add propeller thrust
			//old conditions
			//if (vel_y_old > 0. || angle_old > 0) {		//prototype moving up or prototype pitching down, don't activate thrust
			//prototype still pitching up, don't activate thrust
			//new condition
			//activate after time t
			if (t < activation_time_) {		
				
				propeller_thrust_y = 0.;
				
				propeller_thrust = 0.;
				propeller_thrust_x = 0.;
				
			}
			
			else {
				
				/*if (acc_y_old < 0) {	//provide a thrust opp and equivalent to current decceleration
					
					propeller_thrust_y = -mass_*acc_y_old;
					
				}
				
				else {	//provide a thrust equivalent to current acceleration
					
					propeller_thrust_y = mass_*acc_y_old;
					
				}*/
				//quadratic increase of thrust
				propeller_thrust_y = (t - activation_time_)*(t - activation_time_)*propeller_coeff_;
				
				propeller_thrust = propeller_thrust_y/abs(sin(angle_old));
				
				//can't exceed max_propeller_thrust_
				if (propeller_thrust > max_propeller_thrust_) {
					
					propeller_thrust = max_propeller_thrust_;
					
				}
				
				propeller_thrust_x = propeller_thrust*cos(-angle_old);
				propeller_thrust_y = propeller_thrust*abs(sin(angle_old));
				
			}
			
			//new acc with propeller_thrust_x
			

			x_displ = x_displ_old + delta_x;
			
			//new acc with propeller_thrust_y
			acc_y = accelerationRelaxation_*(f.forceEff()[1]/mass_ + g_[1] + propeller_thrust_y/mass_) + (1. - accelerationRelaxation_)*acc_y_old;
			
			vel_y = vel_y_old + accelerationDamping_*delta_time*(gamma_new_*acc_y + (1. - gamma_new_)*acc_y_old);
			
			delta_y = vel_y*delta_time + accelerationDamping_*delta_time*delta_time*(beta_new_*acc_y + (0.5 - beta_new_)*acc_y_old);

			y_displ = y_displ_old + delta_y;
			
			//Info << "origin old = " << origin_ << " centreOfMass = " << centreOfMass_ << " centreOfRotation = " << motion_.centreOfRotation() << endl;
			
			//origin_ = origin_ + vector(0, delta_y, 0);
			centreOfMass_ = centreOfMass_ + vector(delta_x, delta_y, 0);
			motion_.centreOfRotation() = motion_.centreOfRotation() + vector(delta_x, delta_y, 0);
			
		}
		
		else if (time_.timeIndex() == 1) {
			
			Info << "6dof applied for 1st time" << endl;
			
			//rotational
			angular_acc_z = f.momentEff()[2]/momentOfInertia_[2];
			
			angular_vel_z = angular_vel_z_old + delta_time*angular_acc_z;
			
			angle = angle_old + angular_vel_z_old*delta_time + 0.5*delta_time*delta_time*angular_acc_z;
			
			//x translation
			acc_x = f.forceEff()[0]/mass_ + g_[0];
			
			vel_x = initial_velocity_[0] + delta_time*acc_x;
			
			delta_x = vel_x_old*delta_time + 0.5*delta_time*delta_time*acc_x;
			
			x_displ = x_displ_old + delta_x;
			
			Info << "x_displ = " << x_displ << " vel_x = " << vel_x << " acc_x = " << acc_x << endl;
			
			//y translation
			acc_y = f.forceEff()[1]/mass_ + g_[1];
			
			vel_y = initial_velocity_[1] + delta_time*acc_y;
			
			delta_y = vel_y_old*delta_time + 0.5*delta_time*delta_time*acc_y;
			
			y_displ = y_displ_old + delta_y;
			
			Info << "y_displ = " << y_displ << " vel_y = " << vel_y << " acc_y = " << acc_y << endl;

			//origin_ = origin_ + vector(0, delta_y, 0);
			centreOfMass_ = centreOfMass_ + vector(delta_x, delta_y, 0);
			motion_.centreOfRotation() = motion_.centreOfRotation() + vector(delta_x, delta_y, 0);	
			
		}
		
		if (angleEveryStep_ == true) {
		
			//write into file
			//ofstream MyFile3("current_time_angle_origin_CG.dat");
			//MyFile3 << std::fixed << std::setprecision(8) << t << " " << angle << " " << origin_.component(vector::X) << " " << origin_.component(vector::Y) << " " << origin_.component(vector::Z) << " " << motion_.centreOfRotation().component(vector::X) << " " << motion_.centreOfRotation().component(vector::Y) << " " << motion_.centreOfRotation().component(vector::Z);
			
			//ofstream MyFile3("current_time_angle_thrust_x.dat");
			ofstream MyFile3("current_time_angle_displ_thrust_x.dat");
			MyFile3 << std::fixed << std::setprecision(12) << t << " " << angle << " " << x_displ << " " << y_displ << " " << propeller_thrust_x;

			// Close the file
			MyFile3.close();
			
		}
		
		if (t/fileWriteInterval_ >= saved_no) {
		
			//saving files
			ofstream MyFile("interval_angle_displ_vel_acc_CG.dat");
			MyFile << std::fixed << std::setprecision(8) << t << " " << angle << " " << angular_vel_z << " " << angular_acc_z << " " << x_displ << " " << vel_x << " " << acc_x << " " << y_displ << " " << vel_y << " " << acc_y << " " << motion_.centreOfRotation().component(vector::X) << " " << motion_.centreOfRotation().component(vector::Y) << " " << motion_.centreOfRotation().component(vector::Z);

			// Close the file
			MyFile.close();
						
			saved_no = saved_no + 1;
		}
		//}
		
		
		
		angle_old = angle;
		angular_vel_z_old = angular_vel_z;
		angular_acc_z_old = angular_acc_z;
		
		x_displ_old = x_displ;
		vel_x_old = vel_x;
		acc_x_old = acc_x;
		
		y_displ_old = y_displ;
		vel_y_old = vel_y;
		acc_y_old = acc_y;
		
		displacement[0] = x_displ;
		displacement[1] = y_displ;
		displacement[2] = 0;
		
		saved_time_ = t;
		

		//Info << "current time = " << t << " new angle due to 1dof pitch = " << angle << " new angular velocity = " << angular_vel_z << " new angular acceleration = " << angular_acc_z << endl;
		
		//Info << "6-DoF rigid body motion" << endl;
		Info << "    Pitch angle: " << angle << endl;
		Info << "    Angular velocity: " << angular_vel_z << endl;
		Info << "    Angular acceleration: " << angular_acc_z << endl;
		Info << "    x displacement: " << x_displ << endl;
		Info << "    x velocity: " << vel_x << endl;
		Info << "    x acceleration: " << acc_x << endl;		
		Info << "    y displacement: " << y_displ << endl;
		Info << "    y velocity: " << vel_y << endl;
		Info << "    y acceleration: " << acc_y << endl;
		Info << "    propeller_thrust_x: " << propeller_thrust_x << endl;
		Info << "    propeller_thrust_y: " << propeller_thrust_y << endl;
		Info << "    propeller_thrust: " << propeller_thrust << endl;
		Info << "    origin: " << origin_ << endl;
		Info << "    centreOfMass: " << centreOfMass_ << endl;
		Info << "    centreOfRotation: " << motion_.centreOfRotation() << endl;
		
	}
	
	/*
    quaternion R(axis_, angle);
    //septernion TR(septernion(-origin_)*R*septernion(origin_));
    quaternion R2(1);
	septernion TR(septernion(-displacement)*R2*septernion(-motion_.centreOfRotation())*R*septernion(motion_.centreOfRotation()));    
	*/
	
	//rotation
	quaternion R(axis_, angle);
    septernion TR1(septernion(-origin_)*R*septernion(origin_));quaternion R2(1);
	//then translation
    septernion TR2(septernion(-displacement)*R2);
	septernion TR(TR2*TR1);
    

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::rotatingMotion_3dof_pitch_xy_translate_ext_thrust::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    /*omega_.reset
    (
        Function1<scalar>::New("omega", SBMFCoeffs_, &time_)
    );*/

    return true;
}


// ************************************************************************* //

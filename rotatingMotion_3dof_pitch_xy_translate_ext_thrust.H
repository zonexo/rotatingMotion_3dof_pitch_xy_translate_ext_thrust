/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Class
    Foam::solidBodyMotionFunctions::rotatingMotion_3dof_pitch_xy_translate_ext_thrust

Description
    SolidBodyMotionFvMesh 6DoF motion function.

    The rotation is defined by an origin and axis of rotation and an angular
    speed.

SourceFiles
    rotatingMotion_3dof_pitch_xy_translate_ext_thrust.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingMotion_3dof_pitch_xy_translate_ext_thrust_H
#define rotatingMotion_3dof_pitch_xy_translate_ext_thrust_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"
#include "sixDoFRigidBodyMotion.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rotatingMotion_3dof_pitch_xy_translate_ext_thrust Declaration
\*---------------------------------------------------------------------------*/

class rotatingMotion_3dof_pitch_xy_translate_ext_thrust
:
    public solidBodyMotionFunction
{
		// Private Data
		// 
		dictionary dict_;

        //- Six DoF motion object
        // mutable becuase transformation() function is const.
        mutable sixDoFRigidBodyMotion motion_;

        //- Current time index (used for updating)
		mutable bool firstRun_;
        mutable label startTimeIndex_;	//index to start applying 6dof
		//const label waitToRunIndex_;	//amt of index to wait before apply 6dof 
		//const label waitToRunIndexRamp_;	//amt of index to apply linear bet old and new data to ramp up influence slowly
		const scalar fileWriteInterval_;
		// //- Patches to integrate forces
        // const labelHashSet patchSet_;
		
		mutable scalar angle_old;
		mutable scalar angular_vel_z_old;
		mutable scalar angular_acc_z_old;
		
		//mutable scalar angle_ramp;
		//mutable scalar angular_vel_z_ramp;
		//mutable scalar angular_acc_z_ramp;
		
		wordRes patches_;
		
		//- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;
		
        

        //- Origin of the axis
        const vector origin_;

        //- Axis vector
        const vector axis_;
		const scalar mass_;
		mutable vector centreOfMass_;
		const vector momentOfInertia_;
		const scalar accelerationRelaxation_;
		const scalar accelerationDamping_;
		const vector g_;
		
		
		//- Coefficient for velocity integration (default: 0.5)
        const scalar gamma_;
  
        //- Coefficient for position and orientation integration (default: 0.25)
        const scalar beta_;
         
        //- Linear velocity
		mutable scalar x_displ_old;
		mutable scalar vel_x_old;
		mutable scalar acc_x_old;
		
		mutable scalar y_displ_old;
		mutable scalar vel_y_old;
		mutable scalar acc_y_old;
		
		//mutable scalar y_displ_ramp;
		//mutable scalar vel_y_ramp;
		//mutable scalar acc_y_ramp;
		
		mutable scalar saved_time_;
		mutable label saved_no;
		mutable bool angleEveryStep_;
		const scalar max_propeller_thrust_;
        const scalar activation_pitch_angle_;
        const scalar activation_time_;
        const scalar propeller_coeff_;

        //- Angular velocity (rad/sec)
        //autoPtr<Function1<scalar>> omega_;
        
        


    // Private Member Functions

        //- No copy construct
        rotatingMotion_3dof_pitch_xy_translate_ext_thrust(const rotatingMotion_3dof_pitch_xy_translate_ext_thrust&) = delete;

        //- No copy assignment
        void operator=(const rotatingMotion_3dof_pitch_xy_translate_ext_thrust&) = delete;


public:

    //- Runtime type information
    TypeName("rotatingMotion_3dof_pitch_xy_translate_ext_thrust");


    // Constructors

        //- Construct from components
        rotatingMotion_3dof_pitch_xy_translate_ext_thrust
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new rotatingMotion_3dof_pitch_xy_translate_ext_thrust
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~rotatingMotion_3dof_pitch_xy_translate_ext_thrust() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
		
	
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

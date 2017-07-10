/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    relativeRotor

Description
    Single rotating region

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    relativeRotor.C

\*---------------------------------------------------------------------------*/

#ifndef relativeRotor_H
#define relativeRotor_H

#include "cylindricalCS.H"
#include "Switch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class regionSplit;
class polyTopoChanger;

/*---------------------------------------------------------------------------*\
                        Class relativeRotor Declaration
\*---------------------------------------------------------------------------*/

class relativeRotor
{
    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Coordinate system
        cylindricalCS cs_;

  //const Switch relative_;

  const word relativeRotor_;

  cylindricalCS relativeCs_;

  cylindricalCS& tempCs_;
        // - Rotational speed in rotations per minute (rpm)
  const word motionType_;
  const scalar rpm_;

  scalar startPhase_;
  scalar amplitude_;
  scalar frequency_;

  bool initialize_;

  scalar relativeRpm_;

        //- Name of moving slider patch
        const word movingSliderName_;

        //- Name of static slider patch
        const word staticSliderName_;

        //- Rotating region marker point
        point rotatingRegionMarker_;

        //- Invert motion mask: marked region is stationary and the rest
        //  of the domain is moving
        const Switch invertMotionMask_;

        //- Use topological sliding.  Other option is GGI
        const Switch useTopoSliding_;

        //- Attach-detach action
        Switch attachDetach_;
        
        const word cellSetName_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        relativeRotor(const relativeRotor&);

        //- Disallow default bitwise assignment
        void operator=(const relativeRotor&);


        //- Return coordinate system
        const cylindricalCS& cs() const
        {
            return cs_;
        }

        //- Calculate moving mask
        void calcMovingMask() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;

        //- Clear moving points mask
        void clearPointMask();
  
  Foam::vector getRelativeAngularMotion();

public:

    // Constructors

        //- Construct from  dictionary
        relativeRotor
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~relativeRotor();


    // Member Functions

  void setRelativeMotion
  (
   const coordinateSystem originCs,
   const scalar originRpm
   );

  Foam::word getName() const
  {
    return name_;
  }
  
  Foam::coordinateSystem getCs() const
  {
    return cs_;
  }

  Foam::scalar getRpm() const 
  {
    return rpm_;
  }
  
  Foam::word  getRelativeRotor()
  {
    return relativeRotor_;
  }

  void findAndSetOrigin(PtrList<relativeRotor>& sourceRotors);

        //- Return true if topological sliding is used
        bool useTopoSliding() const
        {
            return useTopoSliding_;
        }

        //- Return true if attach-detach action is used
        bool attachDetach() const
        {
            return attachDetach_;
        }

        //- Return accumulative point motion
        virtual tmp<vectorField> pointMotion() const;

        //- Return accumulative point motion
  Foam::tmp<Foam::vectorField> relativePointMotion();

        //- Add mixer rotor zones
        void addZones
        (
            DynamicList<pointZone*>& pz,
            DynamicList<faceZone*>& fz,
            DynamicList<cellZone*>& cz,
            const regionSplit& rs
        );

        //- Add mixer rotor  modifiers
        void addModifiers
        (
            polyTopoChanger& tc,
            label& nextI
        );

        //- Update topology
        void updateTopology();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

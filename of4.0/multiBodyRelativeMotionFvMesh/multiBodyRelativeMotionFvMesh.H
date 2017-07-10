/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::multiBodyRelativeMotionFvMesh

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    multiBodyRelativeMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef multiBodyRelativeMotionFvMesh_H
#define multiBodyRelativeMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"
#include "solidBodyMotionFunction.H"
//#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiBodyRelativeMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multiBodyRelativeMotionFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        const dictionary dynamicMeshCoeffs_;

  const word refFrame_;

        //- The motion control function
        PtrList<solidBodyMotionFunction> SBMFs_;

  //List<septernion> complexMotionSepternions_;
  labelList relativeZoneIDs_;

        //- The reference points which are transformed
        pointIOField undisplacedPoints_;

        //- Specified cellZones
        labelList zoneIDs_;

  wordList zoneNames_; 

        //- Points to move per cellZone
        labelListList pointIDs_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiBodyRelativeMotionFvMesh(const multiBodyRelativeMotionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const multiBodyRelativeMotionFvMesh&);


public:

    //- Runtime type information
    TypeName("multiBodyRelativeMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        multiBodyRelativeMotionFvMesh(const IOobject& io);


    //- Destructor
    ~multiBodyRelativeMotionFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

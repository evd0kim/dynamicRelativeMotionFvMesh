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
    relativeMotionGgiFvMesh

Description
    Simple mixer mesh using a GGI interface

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    relativeMotionGgiFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef relativeMotionGgiFvMesh_H
#define relativeMotionGgiFvMesh_H

#include "dynamicFvMesh.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class relativeMotionGgiFvMesh Declaration
\*---------------------------------------------------------------------------*/

class relativeMotionGgiFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Coordinate system
        cylindricalCS cs_;

        // - Rotational speed in rotations per minute (rpm)
        scalar rpm_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        relativeMotionGgiFvMesh(const relativeMotionGgiFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const relativeMotionGgiFvMesh&);


        //- Add mixer zones
        void addZonesAndModifiers();

        //- Calculate moving masks
        void calcMovingMasks() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;


public:

    //- Runtime type information
    TypeName("relativeMotionGgiFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit relativeMotionGgiFvMesh(const IOobject& io);


    // Destructor

        virtual ~relativeMotionGgiFvMesh();


    // Member Functions

        //- Return coordinate system
        const cylindricalCS& cs() const
        {
            return cs_;
        }

        //- Update the mesh for both mesh motion
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

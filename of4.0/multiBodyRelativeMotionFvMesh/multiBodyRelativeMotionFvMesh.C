/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "multiBodyRelativeMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "transformField.H"
#include "cellZoneMesh.H"
#include "boolList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiBodyRelativeMotionFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        multiBodyRelativeMotionFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiBodyRelativeMotionFvMesh::multiBodyRelativeMotionFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    undisplacedPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    )
{
    if (undisplacedPoints_.size() != nPoints())
    {
        FatalIOErrorInFunction
        (
            dynamicMeshCoeffs_
        )   << "Read " << undisplacedPoints_.size()
            << " undisplaced points from " << undisplacedPoints_.objectPath()
            << " but the current mesh has " << nPoints()
            << exit(FatalIOError);
    }


    zoneIDs_.setSize(dynamicMeshCoeffs_.size());
    zoneNames_.setSize(dynamicMeshCoeffs_.size());
    SBMFs_.setSize(dynamicMeshCoeffs_.size());
    relativeZoneIDs_.setSize(dynamicMeshCoeffs_.size());
    pointIDs_.setSize(dynamicMeshCoeffs_.size());
    label zoneI = 0;

    forAllConstIter(dictionary, dynamicMeshCoeffs_, iter)
    {
        if (iter().isDict())
        {
            zoneIDs_[zoneI] = cellZones().findZoneID(iter().keyword());
	    zoneNames_[zoneI] = iter().keyword();
	    
            if (zoneIDs_[zoneI] == -1)
            {
                FatalIOErrorInFunction
                (
                    dynamicMeshCoeffs_
                )   << "Cannot find cellZone named " << iter().keyword()
                    << ". Valid zones are " << cellZones().names()
                    << exit(FatalIOError);
            }

            const dictionary& subDict = iter().dict();

	    const word refFrame_
	      (subDict.lookupOrDefault<Foam::word>
	       ("referenceFrame", "none"));

	    if(refFrame_ != "none")
	      {
		relativeZoneIDs_[zoneI] = cellZones().findZoneID(refFrame_);
		if ( relativeZoneIDs_[zoneI] == -1)
		  {
		    FatalIOErrorInFunction
		      (
		       dynamicMeshCoeffs_
		       )   
		      << "Cannot find reference cellZone named " 
		      << refFrame_
		      << ". Valid zones are " 
		      << cellZones().names()
		      << exit(FatalIOError);
		  }
	      }
	    else relativeZoneIDs_[zoneI] = -1;

            SBMFs_.set
            (
                zoneI,
                solidBodyMotionFunction::New(subDict, io.time())
            );

            // Collect points of cell zone.
            const cellZone& cz = cellZones()[zoneIDs_[zoneI]];

            boolList movePts(nPoints(), false);

            forAll(cz, i)
            {
                label celli = cz[i];
                const cell& c = cells()[celli];
                forAll(c, j)
                {
                    const face& f = faces()[c[j]];
                    forAll(f, k)
                    {
                        label pointi = f[k];
                        movePts[pointi] = true;
                    }
                }
            }

            syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

            DynamicList<label> ptIDs(nPoints());
            forAll(movePts, i)
            {
                if (movePts[i])
                {
                    ptIDs.append(i);
                }
            }

            pointIDs_[zoneI].transfer(ptIDs);

            Info<< "Applying solid body motion " 
		<< SBMFs_[zoneI].type()
                << " to " 
		<< pointIDs_[zoneI].size() << " points of cellZone "
                << iter().keyword() << endl;
	    
	    if(relativeZoneIDs_[zoneI] != -1)
	      {
		Info<< "The motion CS has relative reference frame " 
		    << zoneNames_[relativeZoneIDs_[zoneI]]<< endl;
	      }
	    else
	      {
		Info<< "The motion CS has global reference frame " 
		    << endl;
	      }
            zoneI++;
        }
    }
    zoneIDs_.setSize(zoneI);
    zoneNames_.setSize(zoneI);
    SBMFs_.setSize(zoneI);
    relativeZoneIDs_.setSize(zoneI);
    pointIDs_.setSize(zoneI);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiBodyRelativeMotionFvMesh::~multiBodyRelativeMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiBodyRelativeMotionFvMesh::update()
{
    static bool hasWarned = false;

    pointField transformedPts(undisplacedPoints_);

    forAll(zoneIDs_, i)
    {
        const labelList& zonePoints = pointIDs_[i];
	
	//initialization
	vector zV(0,0,0);
	quaternion R(1);
	septernion TR(septernion(zV)*R);
	
	//first, reference
	if(relativeZoneIDs_[i]!=-1)
	  {
	    label idx = findIndex(zoneIDs_, relativeZoneIDs_[i]);
	    TR *= 
	      SBMFs_[idx].transformation();
	  }
	//second, relative
	TR *= SBMFs_[i].transformation();

        UIndirectList<point>(transformedPts, zonePoints) =
            transformPoints
            (
	     TR,
	     pointField(transformedPts, zonePoints)
            );
    }

    fvMesh::movePoints(transformedPts);

    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningInFunction
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //

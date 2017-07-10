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

\*---------------------------------------------------------------------------*/

#include "relativeRotor.H"
#include "regionSplit.H"
#include "polyTopoChanger.H"
#include "slidingInterface.H"
#include "Time.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::relativeRotor::addZones
(
    DynamicList<pointZone*>& pz,
    DynamicList<faceZone*>& fz,
    DynamicList<cellZone*>& cz,
    const regionSplit& rs
)
{
    // Get the region of the cell containing the origin.
    label originRegion = rs[mesh_.findNearestCell(rotatingRegionMarker_)];

    labelList movingCells(mesh_.nCells());
    label nMovingCells = 0;

    forAll(rs, cellI)
    {
        if (rs[cellI] == originRegion)
        {
            movingCells[nMovingCells] = cellI;
            nMovingCells++;
        }
    }

    movingCells.setSize(nMovingCells);
    Info<< "Number of moving cells for " << name_ << ": "
        << nMovingCells << endl;

    cz.append
    (
        new cellZone
        (
            "movingCellsZone" + name_,
            movingCells,
            cz.size(),
            mesh_.cellZones()
        )
    );

    if (useTopoSliding_)
    {
        // Add an empty zone for cut points
        pz.append
        (
            new pointZone
            (
                "cutPointZone" + name_,
                labelList(0),
                pz.size(),
                mesh_.pointZones()
            )
        );


        // Do face zones for slider

        // Moving slider
        label movingSliderIndex =
            mesh_.boundaryMesh().findPatchID(movingSliderName_);

        if (movingSliderIndex < 0)
        {
            FatalErrorIn("void relativeRotor::addZones(...) const")
                << "Moving slider patch not found in boundary"
                << abort(FatalError);
        }

        label staticSliderIndex =
            mesh_.boundaryMesh().findPatchID(staticSliderName_);

        if (staticSliderIndex < 0)
        {
            FatalErrorIn("void relativeRotor::addZones(...) const")
                << "Static slider patch not found in boundary"
                << abort(FatalError);
        }

        const polyPatch& movingSlider =
            mesh_.boundaryMesh()[movingSliderIndex];

        labelList isf(movingSlider.size());

        forAll (isf, i)
        {
            isf[i] = movingSlider.start() + i;
        }

        fz.append
        (
            new faceZone
            (
                movingSliderName_ + "Zone" + name_,
                isf,
                boolList(movingSlider.size(), false),
                fz.size(),
                mesh_.faceZones()
            )
        );

        // Static slider
        const polyPatch& staticSlider =
            mesh_.boundaryMesh()[staticSliderIndex];

        labelList osf(staticSlider.size());

        forAll (osf, i)
        {
            osf[i] = staticSlider.start() + i;
        }

        fz.append
        (
            new faceZone
            (
                staticSliderName_ + "Zone" + name_,
                osf,
                boolList(staticSlider.size(), false),
                fz.size(),
                mesh_.faceZones()
            )
        );

        // Add empty zone for cut faces
        fz.append
        (
            new faceZone
            (
                "cutFaceZone" + name_,
                labelList(0),
                boolList(0, false),
                fz.size(),
                mesh_.faceZones()
            )
        );
    }
}


void Foam::relativeRotor::addModifiers
(
    polyTopoChanger& tc,
    label& nextI
)
{
    // Add a topology modifier
    if (useTopoSliding())
    {
        Info << "Adding topology modifier for rotor " << name_ << endl;

        tc.set
        (
            nextI,
            new slidingInterface
            (
                "mixerSlider"  + name_,
                nextI,
                tc,
                staticSliderName_ + "Zone"  + name_,
                movingSliderName_ + "Zone" + name_,
                "cutPointZone" + name_,
                "cutFaceZone" + name_,
                staticSliderName_,
                movingSliderName_,
                slidingInterface::INTEGRAL,   // Edge matching algorithm
                attachDetach_,                // Attach-detach action
                intersection::VISIBLE         // Projection algorithm
            )
        );

        nextI++;
    }
}


void Foam::relativeRotor::calcMovingMask() const
{
    if (movingPointsMaskPtr_)
    {
        FatalErrorIn("void relativeRotor::calcMovingMask() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = 
      new scalarField(mesh_.allPoints().size(), 0.);

    scalarField& movingPointsMask = *movingPointsMaskPtr_;
    
    Info()<<"Attempt to check field size"<<endl;
    Info()<<"Mask size: "<<movingPointsMask.size()<<endl;

    const cellList& c = mesh_.cells();
    const faceList& f = mesh_.allFaces();

    //problem is here
    if(cellSetName_=="none")
    {
        FatalErrorIn("void relativeRotor::calcMovingMask() const")
	  << "cellZone should be set for the current implementation"
	  << "of the relative motion"
	  << abort(FatalError);
    }

    const labelList& cellAddr = mesh_.cellZones()
        [mesh_.cellZones().findZoneID(cellSetName_)];// "movingCellsZone" + name_)]; // 
    
    //const labelList& cellAddr =
    //    cellZones()[cellZones().findZoneID("movingCells")];

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1.;
            }
        }
    }
    // Attempt to enforce motion on sliders if zones exist
    const label msI =
        mesh_.faceZones().findZoneID(movingSliderName_ + "Zone" + name_);

    if (msI > -1)
    {
        const labelList& movingSliderAddr = mesh_.faceZones()[msI];

        forAll (movingSliderAddr, faceI)
        {
            const face& curFace = f[movingSliderAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }
    }

    const label ssI =
        mesh_.faceZones().findZoneID(staticSliderName_ + "Zone" + name_);

    if (ssI > -1)
    {
        const labelList& staticSliderAddr = mesh_.faceZones()[ssI];

        forAll (staticSliderAddr, faceI)
        {
            const face& curFace = f[staticSliderAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 0;
            }
        }
    }
}


// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::relativeRotor::movingPointsMask() const
{
    if (!movingPointsMaskPtr_)
    {
        calcMovingMask();
    }

    return *movingPointsMaskPtr_;
}


void Foam::relativeRotor::clearPointMask()
{
    deleteDemandDrivenData(movingPointsMaskPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeRotor::relativeRotor
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    cs_
    (
        "coordinateSystem",
        dict.subDict("coordinateSystem")
    ),
    //relative_(dict.lookupOrDefault<bool>("relative", false)),
    relativeRotor_
    (
     dict.lookupOrDefault<word>("relativeRotor", "none")
     ),
    relativeCs_(cs_),
    tempCs_(cs_),
    motionType_(dict.lookupOrDefault<word>("motionType", "rotating")),
    rpm_(readScalar(dict.lookup("rpm"))),
    startPhase_(0.0),
    amplitude_(0.0),
    frequency_(0.0),
    relativeRpm_(rpm_),
    movingSliderName_(dict.lookup("movingPatch")),
    staticSliderName_(dict.lookup("staticPatch")),
    rotatingRegionMarker_
    (
        dict.lookupOrDefault<point>("rotatingRegionMarker", cs_.origin())
    ),
    invertMotionMask_
    (
        dict.lookupOrDefault<bool>("invertMotionMask", false)
    ),
    useTopoSliding_(dict.lookup("useTopoSliding")),
    attachDetach_(dict.lookupOrDefault<bool>("attachDetach", true)),
    cellSetName_(dict.lookupOrDefault<word>("cellSetName", "none")),
    movingPointsMaskPtr_(NULL)
{
  if(motionType_ == "oscillating")
    {
      dictionary coeffsTemp = dict.subDict("oscillatingMotionCoeffs");
      amplitude_ = readScalar(coeffsTemp.lookup("amplitude"));
      startPhase_ = readScalar(coeffsTemp.lookup("phase"));
      frequency_ = readScalar(coeffsTemp.lookup("frequency"));
      initialize_ = readBool(coeffsTemp.lookup("initPhase"));
    }
    // Make sure the coordinate system does not operate in degrees
    // Bug fix, HJ, 3/Oct/2011
    if (!cs_.inDegrees())
    {
        WarningIn
        (
            "relativeRotor::relativeRotor\n"
            "(\n"
            "    const word& name,\n"
            "    const polyMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Mixer coordinate system is set to operate in radians.  "
            << "Changing to rad for correct calculation of angular velocity."
            << nl
            << "To remove this message please add entry" << nl << nl
            << "inDegrees true;" << nl << nl
            << "to the specification of the coordinate system"
            << endl;

        cs_.inDegrees() = true;
    }

    Info<< "Rotor " << name << ":" << nl
        << "    origin      : " << cs().origin() << nl
        << "    axis        : " << cs().axis() << nl
        << "    rpm         : " << rpm_ << nl
        << "    invert mask : " << invertMotionMask_ << nl
        << "    topo sliding: " << useTopoSliding_ << endl;

    if (useTopoSliding_)
    {
        Info<< "    attach-detach: " << attachDetach_ << endl;
    }

    /*if (relativeRotor_ == "none")
    {
        FatalErrorIn("relativeRotor::relativeRotor()")
            << "There is no relative rotor provided in the dictionary"
            << abort(FatalError);
    }
    */
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeRotor::~relativeRotor()
{
    clearPointMask();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::relativeRotor::pointMotion() const
{
    // Rotational speed needs to be converted from rpm
    scalarField mpm = movingPointsMask();

    if (invertMotionMask_)
    {
        Info << "Inverting motion mask" << endl;
        mpm = 1 - mpm;
    }

    vector globalTangentialVector
      (0, 
       relativeRpm_*360.0*mesh_.time().deltaT().value()/60.0
       ,0);

    tmp<Foam::vectorField> globalTransformVector
      (
       relativeCs_.globalPosition
       (
        // Motion vector in cylindrical coordinate system (x theta z)
        relativeCs_.localPosition(mesh_.allPoints())
	+ globalTangentialVector*mpm
	)
	- mesh_.allPoints()
       );

    return globalTransformVector;
}

Foam::tmp<Foam::vectorField> Foam::relativeRotor::relativePointMotion()
{
  if( relativeRotor_!= "none")
    {	
      
      // Rotational speed needs to be converted from rpm
      scalarField mpm = movingPointsMask();
      
      if (invertMotionMask_)
	{
	  Info << "Inverting motion mask" << endl;
	  mpm = 1 - mpm;
	}
      
      vector globalTangentialVector
	(0, 
	 relativeRpm_*360.0*mesh_.time().deltaT().value()/60.0,
	 0);
           
      if(true)
	{
	  Info()<<"Debug: "<<name_
		<<" local coordinate system origin "
		<<tempCs_.origin()<<nl
		<<"\tDirection vector"<<tempCs_.direction()<<nl
		<<"\tWhole cs info "<<tempCs_<<endl;
	}

      point tOrigin = //cs_.origin();
	relativeCs_.globalPosition
	(
	 relativeCs_.localPosition(tempCs_.origin())
	 + globalTangentialVector
	 ); 
	
      vector tAxis = tempCs_.axis();
      vector tDirn = Foam::vector(tOrigin);

      cylindricalCS relativeCs("cylindrical", 
			   tOrigin, 
			   tAxis, 
			   tDirn, 
			   tempCs_.inDegrees());
      if(true)
	{
	  Info()<<"Debug: transformed coordinate system origin "
		<<tOrigin<<nl
		<<"\tDirection vector"<<tDirn<<nl
		<<"\tWhole cs info"<<relativeCs<<endl;
	}

      Foam::vector localTangentialVector = getRelativeAngularMotion();
     
      tmp<Foam::vectorField> localTransformVector 
	(
	 relativeCs.globalPosition
	 (
	  relativeCs.localPosition(mesh_.allPoints())
	  + localTangentialVector*mpm
	  )
	 - mesh_.allPoints());

      tempCs_ = relativeCs;
      
      return localTransformVector; 
    }
  
  return tmp<Foam::vectorField>
    (Foam::vectorField(mesh_.allPoints().size(), vector::zero));

}

Foam::vector Foam::relativeRotor::getRelativeAngularMotion()  
{
  Foam::vector 
    returnVector(.0, 
		 rpm_*360.0*mesh_.time().deltaT().value()/60.0,
		 .0);
  
  if(motionType_ == "oscillating")
    {
      //const scalar deg2Rad = mathematicalConstant::pi/180.0;
      Foam::scalar phase =  
	frequency_*mesh_.time().value()+startPhase_;
      //phase *= deg2Rad;
      returnVector[1] = 
	amplitude_*frequency_*cos(phase)*
	mesh_.time().deltaT().value();

      if(initialize_)
      {
	returnVector[1] += amplitude_*sin(startPhase_);
	initialize_ = false;
      };
      /*
      Info()<<"DEBUG:"<<nl
	    <<"freq="<<frequency_<<nl
	    <<"time="<<mesh_.time().value()<<nl
	    <<"phase0="<<startPhase_<<nl
	    <<"returnVector="<<returnVector[1]<<endl;
      */
    }
  return returnVector;
}

void Foam::relativeRotor::updateTopology()
{
    clearPointMask();
}

void Foam::relativeRotor::setRelativeMotion(const coordinateSystem originCs, const scalar originRpm)
{
  relativeCs_ = originCs;
  relativeRpm_ = originRpm;
}

void Foam::relativeRotor::findAndSetOrigin(PtrList<relativeRotor>& sourceRotors)
{
  forAll(sourceRotors, rotorI)
    {
      if(sourceRotors[rotorI].getName() == relativeRotor_)
	{
	  Info()<<"\tDebug: Relative rotor has been found"<<endl;
	  Info()<<"\tOrigin name "<<sourceRotors[rotorI].getName()<<nl
		<<"\tCurrent rotor name "<<name_<<endl;
	  relativeCs_ = sourceRotors[rotorI].getCs();
	  relativeRpm_ = sourceRotors[rotorI].getRpm();   
	}
    }
}

// ************************************************************************* //

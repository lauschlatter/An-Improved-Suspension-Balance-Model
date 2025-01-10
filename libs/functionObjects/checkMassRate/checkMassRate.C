/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    LTFD.
    Laboratory of Thermo-Fluid Dynamics (LTFD)
    Chemical Engineering Program (PEQ)
    Alberto Luiz Coimbra Institute for Graduate Studies and Research in Engineering (COPPE)
    Federal University of Rio de Janeiro Rio de Janeiro, Brazil
    
Written by:
    Lauren Schlatter Fernandes
    laurenfernandes@peq.coppe.ufrj.br

    Gabriel Gonçalves da Silva Ferreira
    gferreira@peq.coppe.ufrj.br


\*----------------------------------------------------------------------------*/

#include "checkMassRate.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "fvCFD.H"
#include "suspensionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(checkMassRate, 0);

        addToRunTimeSelectionTable
        (
            functionObject,
            checkMassRate,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::checkMassRate::checkMassRate
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    inletPatchID_(),
    outletPatchID_(),
    historyFilePtr_(NULL)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    word inletPatchName;
    word outletPatchName;
    dict.readIfPresent("inletPatch", inletPatchName);
    dict.readIfPresent("outletPatch", outletPatchName);
    inletPatchID_ = mesh.boundaryMesh().findPatchID(inletPatchName);
    outletPatchID_ = mesh.boundaryMesh().findPatchID(outletPatchName);

    if ((inletPatchID_ == -1) || (outletPatchID_ == -1))
    {
        FatalErrorIn("Foam::checkMassRate::read()")
                    << "inletPatch and outletPatch names must be"
                    << " correctly informed" << exit(FatalError);
    }

    if (Pstream::master())
    {
        if (!time_.processorCase())
        {
            mkDir
            (
                time_.path()
                /"history"
                /time_.timeName()
            );

            historyFilePtr_ =
                new OFstream
                (
                    time_.path()
                    /"history"
                    /time_.timeName()
                    /"checkMassRate.dat"
                );
        }
        else
        {
            mkDir
            (
                time_.path()/".."/"history"
                /time_.timeName()
            );

            historyFilePtr_ =
                new OFstream
                (
                    time_.path()/".."
                    /"history"
                    /time_.timeName()
                    /"checkMassRate.dat"
                );
        }

        (*historyFilePtr_)
            << "Time" << tab 
            << "inFlux" << tab
            << "outFlux" << tab
            << "Diff" << tab
            << "pIn" << tab
            << "pOut" << tab
            << "delta p" << tab
            << "End" << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::checkMassRate::start()
{
    Info<< "-----------------------------------------" << endl;
    Info<< "functionObject::checkMassRate::start()" << endl;
    Info<< "inletPatchID = " << inletPatchID_
        << "\toutletPatchID = " << outletPatchID_ << endl;
    Info<< "-----------------------------------------" << endl;

    return true;
}


bool Foam::functionObjects::checkMassRate::execute()
{
    Info<< "-----------------------------------------" << endl;
    Info<< "functionObject::checkMassRate::execute()" << endl;
    Info<< "-----------------------------------------" << endl;
    
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    const suspensionModel& mixture =
        mesh.lookupObject<suspensionModel>("suspensionProperties");

// Total mass flux in multiphase system at the domain's entrance and exit
    scalar influx(0);
    scalar outflux(0);
    scalar pinave(0);
    scalar poutave(0);
    scalar ain(0);
    scalar aout(0);

    const scalarField& cIn = mixture.c().boundaryField()[inletPatchID_];
    const scalarField& cOut = mixture.c().boundaryField()[outletPatchID_];
    const scalarField& pIn = mixture.p().boundaryField()[inletPatchID_];
    const scalarField& pOut = mixture.p().boundaryField()[outletPatchID_];
    const scalarField& phiIn = mixture.phi().boundaryField()[inletPatchID_];
    const scalarField& phiOut = mixture.phi().boundaryField()[outletPatchID_];
    const scalarField& aIn = mesh.magSf().boundaryField()[inletPatchID_];
    const scalarField& aOut = mesh.magSf().boundaryField()[outletPatchID_];

    forAll(cIn, faceI)
    {
        influx += phiIn[faceI]*cIn[faceI];
        pinave += pIn[faceI]*aIn[faceI];
        ain += aIn[faceI];
    }

    forAll(cOut, faceI)
    {
        outflux += phiOut[faceI]*cOut[faceI];
        poutave += pOut[faceI]*aOut[faceI];
        aout += aOut[faceI];
    }

    reduce(influx, sumOp<scalar>());
    reduce(outflux, sumOp<scalar>());

    reduce(pinave, sumOp<scalar>());
    reduce(poutave, sumOp<scalar>());

    reduce(ain, sumOp<scalar>());
    reduce(aout, sumOp<scalar>());

    pinave = pinave/ain;
    poutave = poutave/aout;

    Info<< "Particle's Influx = " << influx << " m^3 / s"
        << "\tParticle's Outflux = " << outflux << " m^3 / s"
        << "\t Diff = " << influx + outflux << " m^3 / s" << endl;
    Info<< "p in ave = " << pinave
        << "\tp out ave = " << poutave
        << "\t -delta p = " << pinave - poutave << endl;
    Info<< "----------------------------------------------------------"
        << endl
        << endl;

    if (Pstream::master())
    {
        (*historyFilePtr_)
            << time_.value() << tab 
            << influx << tab
            << outflux << tab
            << influx + outflux << tab
            << pinave << tab
            << poutave << tab
            << pinave - poutave << tab;

        (*historyFilePtr_)
            << endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::checkMassRate::write()
{
    return true;
}


bool Foam::functionObjects::checkMassRate::read(const dictionary& dict)
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    word inletPatchName;
    word outletPatchName;

    dict.readIfPresent("inletPatch", inletPatchName);
    dict.readIfPresent("outletPatch", outletPatchName);

    inletPatchID_ = mesh.boundaryMesh().findPatchID(inletPatchName);
    outletPatchID_ = mesh.boundaryMesh().findPatchID(outletPatchName);

    if ((inletPatchID_ == -1) || (outletPatchID_ == -1))
    {
        FatalErrorIn("Foam::checkMassRate::read()")
                    << "inletPatch and outletPatch names must be correctly informed"
                    << exit(FatalError);
    }

    return true;
}

// ************************************************************************* //

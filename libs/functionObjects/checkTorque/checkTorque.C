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

#include "checkTorque.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(checkTorque, 0);
    addToRunTimeSelectionTable(functionObject, checkTorque, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::checkTorque::checkTorque
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
    CofR_(vector::zero),
    patchSet_(),
    historyFilePtr_(NULL)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

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
                    /"checkTorque.dat"
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
                    /"checkTorque.dat"
                );
        }

        (*historyFilePtr_)
            << "Time" << tab; 

        (*historyFilePtr_)
            << "Torque mix.x" << tab
            << "Torque mix.y" << tab
            << "Torque mix.z";

        (*historyFilePtr_)
            << "End" << endl;
    }

    CofR_ = dict.lookup("CofR");
    patchSet_ =
        mesh.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));    

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::checkTorque::start()
{
    return true;
}


bool Foam::functionObjects::checkTorque::execute()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    if (Pstream::master())
        {
            (*historyFilePtr_)
                << time_.value() << tab;
        }

// STEP 1 : Calculate the stress tensor to be evaluated
    const volTensorField& Tensor = mesh.lookupObject<volTensorField>("Sigma");

// STEP 2: Calculate the torque
    vector TangentMTensor(vector::zero);

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        vectorField Md = mesh.C().boundaryField()[patchi] - CofR_;
        scalarField sA = mag(mesh.Sf().boundaryField()[patchi]);

        vectorField TensordA =
            (mesh.Sf().boundaryField()[patchi] & Tensor.boundaryField()[patchi])
            /sA;

        vectorField fNtensor =
            (mesh.Sf().boundaryField()[patchi]/sA)
            *(mesh.Sf().boundaryField()[patchi] & TensordA);

        vectorField fTtensor = sA*TensordA - fNtensor;

        TangentMTensor += sum(Md ^ fTtensor);
        TangentMTensor = -1*TangentMTensor;
    }

// Parallelize 'sum' operator
    reduce(TangentMTensor,sumOp<vector>());

// STEP 3: Write the results in 'history' file
    if (Pstream::master())
    {
        (*historyFilePtr_)
            << TangentMTensor.x() << tab
            << TangentMTensor.y() << tab
            << TangentMTensor.z() << endl;            
    }

    return true;
}


bool Foam::functionObjects::checkTorque::write()
{
    return true;
}


bool Foam::functionObjects::checkTorque::read(const dictionary& dict)
{
    Info<< "-----------------------------------------" << endl;
    Info<< "functionObject::checkTorque::read()" << endl;
    Info<< "-----------------------------------------" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //

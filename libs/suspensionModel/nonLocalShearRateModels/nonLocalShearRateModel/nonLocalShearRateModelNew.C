/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "nonLocalShearRateModel.H"
#include "volFields.H"

// ************************************************************************* //

Foam::autoPtr<Foam::nonLocalShearRateModel> Foam::nonLocalShearRateModel::New
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
{
    const word modelType(dict.lookup("nonLocalShearRateModel"));

    Info<< "Selecting non local shear rate model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown nonLocalShearRateModel type "
            << modelType << nl << nl
            << "Valid nonLocalShearRateModel are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<nonLocalShearRateModel>
        (cstrIter()(mesh, U, c, dict));
}


// ************************************************************************* //

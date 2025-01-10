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

Class
    Foam::kGradAlpha

SourceFiles
    kGradAlpha.C

Description
    Foam::kGradAlpha

\*---------------------------------------------------------------------------*/

#ifndef kGradAlpha_H
#define kGradAlpha_H

#include "volFields.H"
#include "solidsPressureModel.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidsPressureModels
{
/*---------------------------------------------------------------------------*\
                          Class kGradAlpha Declaration
\*---------------------------------------------------------------------------*/

class kGradAlpha : public solidsPressureModel
{

    protected:

    //- Protected data
    dictionary kGradAlphaDict_;
    dimensionedScalar k_;


    //- Disallow copy construct
    kGradAlpha(const kGradAlpha&);


    public:

    //- Runtime type information
    TypeName("kGradAlpha");


    //- Construct from components
    kGradAlpha
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~kGradAlpha()
    {}


    //- Member Functions
    virtual void updateSolidsPressure();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace solidsPressureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

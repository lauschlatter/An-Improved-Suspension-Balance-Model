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
    Foam::constantNonLocalShearRate

SourceFiles
    constantNonLocalShearRate.C
    
Description
    Foam::constantNonLocalShearRate

\*---------------------------------------------------------------------------*/

#ifndef constantNonLocalShearRate_H
#define constantNonLocalShearRate_H

#include "volFields.H"
#include "nonLocalShearRateModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nonLocalShearRateModels
{
/*---------------------------------------------------------------------------*\
                Class constantNonLocalShearRate Declaration
\*---------------------------------------------------------------------------*/

class constantNonLocalShearRate : public nonLocalShearRateModel
{

    protected:

    //- Protected data
    dictionary constantNonLocalShearRateDict_;
    dimensionedScalar epsilon_;


    //- Disallow copy construct
    constantNonLocalShearRate(const constantNonLocalShearRate&);


    public:

    //- Runtime type information
    TypeName("constantNonLocalShearRate");


    //- Construct from components
    constantNonLocalShearRate
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~constantNonLocalShearRate()
    {}


    //- Member Functions
    virtual void updateNonLocalShearRate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace nonLocalShearRateModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

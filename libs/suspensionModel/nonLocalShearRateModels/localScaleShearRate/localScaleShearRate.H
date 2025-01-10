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
    Foam::localScaleShearRate

SourceFiles
    localScaleShearRate.C
    
Description
    Foam::localScaleShearRate

\*---------------------------------------------------------------------------*/

#ifndef localScaleShearRate_H
#define localScaleShearRate_H

#include "volFields.H"
#include "nonLocalShearRateModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nonLocalShearRateModels
{
/*---------------------------------------------------------------------------*\
                      Class localScaleShearRate Declaration
\*---------------------------------------------------------------------------*/

class localScaleShearRate : public nonLocalShearRateModel
{

    protected:

    //- Protected data
    dictionary localScaleShearRateDict_;
    dimensionedScalar a_;
    dimensionedScalar k_;


    //- Disallow copy construct
    localScaleShearRate(const localScaleShearRate&);


    public:

    //- Runtime type information
    TypeName("localScaleShearRate");


    //- Construct from components
    localScaleShearRate
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~localScaleShearRate()
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

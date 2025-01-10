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
    Foam::expGradAlpha

SourceFiles
    expGradAlpha.C
    
Description
    Foam::expGradAlpha

\*---------------------------------------------------------------------------*/

#ifndef expGradAlpha_H
#define expGradAlpha_H

#include "volFields.H"
#include "solidsPressureModel.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidsPressureModels
{
/*---------------------------------------------------------------------------*\
                          Class expGradAlpha Declaration
\*---------------------------------------------------------------------------*/

class expGradAlpha : public solidsPressureModel
{

    protected:

    //- Protected data
    dictionary expGradAlphaDict_;
    dimensionedScalar preAlphaExp_;
    dimensionedScalar expMax_;
    dimensionedScalar g0_;
    dimensionedScalar cMax_;


    //- Disallow copy construct
    expGradAlpha(const expGradAlpha&);


    public:

    //- Runtime type information
    TypeName("expGradAlpha");


    //- Construct from components
    expGradAlpha
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~expGradAlpha()
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

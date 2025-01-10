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
    Foam::solidsPressureModel

Description
    solidsPressureModel module

\*---------------------------------------------------------------------------*/

#ifndef solidsPressureModel_H
#define solidsPressureModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class solidsPressureModel Declaration
\*---------------------------------------------------------------------------*/

class solidsPressureModel
{

    protected:

    //- Protected data
    const fvMesh& mesh_;
    const volVectorField& U_;
    const volScalarField& c_;
    volVectorField Fp_;
    surfaceScalarField Fpf_;


    //- Disallow copy construct
    solidsPressureModel(const solidsPressureModel&);


    public:

    //- Runtime type information
    TypeName("solidsPressureModel");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidsPressureModel,
        dictionary,
        (
            const fvMesh& mesh,
            const volVectorField& U,
            const volScalarField& c,
            const dictionary& dict
        ),
        (mesh, U, c, dict)
    );


    //- Return a reference to the selected viscosity model
    static autoPtr<solidsPressureModel> New
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Construct from components
    solidsPressureModel
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~solidsPressureModel()
    {}


    //- Member Functions
    virtual void updateSolidsPressure() = 0;
    const volVectorField& Fp() const
    {
        return Fp_;
    }
    volVectorField& Fp()
    {
        return Fp_;
    }
    const surfaceScalarField& Fpf() const
    {
        return Fpf_;
    }
    surfaceScalarField& Fpf()
    {
        return Fpf_;
    }

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

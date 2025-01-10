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
    Foam::shearViscosityModel

Description
    shearViscosityModel module

\*---------------------------------------------------------------------------*/

#ifndef shearViscosityModel_H
#define shearViscosityModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class shearViscosityModel Declaration
\*---------------------------------------------------------------------------*/

class shearViscosityModel
{

    protected:

    //- Protected data
    const fvMesh& mesh_;
    const volVectorField& U_;
    const volScalarField& c_;
    dimensionedScalar cMax_;
    dimensionedScalar nu0_;
    volScalarField nus_;
    volScalarField nup_;


    //- Disallow copy construct
    shearViscosityModel(const shearViscosityModel&);


    public:

    //- Runtime type information
    TypeName("shearViscosityModel");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        shearViscosityModel,
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
    static autoPtr<shearViscosityModel> New
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Construct from components
    shearViscosityModel
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~shearViscosityModel()
    {}


    //- Member Functions
    virtual void updateShearViscosity() = 0;
    
    const dimensionedScalar& nu0() const
    {
        return nu0_;
    }
    dimensionedScalar& nu0()
    {
        return nu0_;
    }
    const dimensionedScalar& cMax() const
    {
        return cMax_;
    }
    dimensionedScalar& cMax()
    {
        return cMax_;
    }
    const volScalarField& nus() const
    {
        return nus_;
    }
    volScalarField& nus()
    {
        return nus_;
    }
    const volScalarField& nup() const
    {
        return nup_;
    }
    volScalarField& nup()
    {
        return nup_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

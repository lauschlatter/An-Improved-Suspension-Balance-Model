/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2015 ESI-OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    

\*---------------------------------------------------------------------------*/

#include "suspensionModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
suspensionModel::suspensionModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& c,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "suspensionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    suspensionDict_(dict.subDict("suspensionProperties")),
    p_(p),
    U_(U),
    a_(dict.subDict("suspensionProperties").lookup("a") ),
    nu0_(dict.subDict("suspensionProperties").lookup("nu0")),
    cMax_(dict.subDict("suspensionProperties").lookup("cMax")),
    J_
    (
        IOobject
        (
            "J",
            U.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Sigma_
    (
        IOobject
        (
            "Sigma",
            U.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "Sigma", dimensionSet(0,2,-2,0,0,0,0), tensor(0,0,0,0,0,0,0,0,0)
        ),
        "calculated"
    ),
    gamma_
    (
        IOobject
        (
            "gamma",
            U.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("gamma", dimensionSet(0,0,-1,0,0,0,0), 0),
        "calculated"
    ),
    gammaNL_
    (
        IOobject
        (
            "gammaNL",
            U.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("gammaNL", dimensionSet(0,0,-1,0,0,0,0), 0),
        "calculated"
    ),
    c_(c),
    phi_(phi),
    I_
    (
        "I",
        dimensionSet(0,0,0,0,0,0,0),
        tensor(1,0,0,0,1,0,0,0,1)
    ),
    divSigmaL1_
    (
        IOobject
        (
            "divSigmaL1",
            U.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE//IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "divSigmaL1", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0)
        ),
        "calculated"
    ),
    divSigmaL1f_
    (
        IOobject
        (
            "divSigmaL1f",
            U.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE//IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("divSigmaL1f", dimensionSet(0,3,-2,0,0,0,0), 0),
        "calculated"
    ),
    divSigmaL2f_
    (
        IOobject
        (
            "divSigmaL2f",
            U.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE//IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("divSigmaL2f", dimensionSet(0,3,-2,0,0,0,0), 0),
        "calculated"
    ),
    nonLocalShearRateModelPtr_
    (
        nonLocalShearRateModel::New
        (
            mesh, U, c, dict.subDict("suspensionProperties")
        )
    ),
    anisotropicStressTensorModelPtr_
    (
        anisotropicStressTensorModel::New
        (
            mesh, U, c, dict.subDict("suspensionProperties")
        )
    ),
    shearViscosityModelPtr_
    (
        shearViscosityModel::New
        (
            mesh, U, c, dict.subDict("suspensionProperties")
        )
    ),
    normalViscosityModelPtr_
    (
        normalViscosityModel::New
        (
            mesh, U, c, dict.subDict("suspensionProperties")
        )
    ),
    hindranceModelPtr_
    (
        hindranceModel::New
        (
            mesh, U, c, dict.subDict("suspensionProperties")
        )
    ),
    solidsPressureModelPtr_
    (
        solidsPressureModel::New
        (
            mesh, U, c, dict.subDict("suspensionProperties")
        )
    )
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //
        
    void Foam::suspensionModel::update()
    {
        Info<< "Updating suspension properties" << endl;

        hindranceModelPtr_->updateHindrance();
        shearViscosityModelPtr_->updateShearViscosity();
        normalViscosityModelPtr_->updateNormalViscosity();

        nonLocalShearRateModelPtr_->updateNonLocalShearRate();        
        anisotropicStressTensorModelPtr_->updateAnisotropicStressTensor();        
        solidsPressureModelPtr_->updateSolidsPressure();

        Sigma_ = SigmaP() + SigmaF();
        updateGamma();
    }


    void Foam::suspensionModel::updateGamma()
    {
        gamma_ = sqrt((2.0*(E() && E())));
        gammaNL_ = nonLocalShearRateModelPtr_->gammaNL();
    }


    void Foam::suspensionModel::updateDivSigmaL1f()
    {
        divSigmaL1_ =
            -dNunDc()*(gammaDot()*(fvc::grad(c_) & Q())
            + gammaNL()*fvc::grad(c_) );

        divSigmaL1f_ =
        (
            -fvc::interpolate(gammaDot()*dNunDc())*fvc::snGrad(c_)*Qf()
            - fvc::interpolate(gammaNL()*dNunDc())*fvc::snGrad(c_)*mesh_.magSf()
        );
    }


    void Foam::suspensionModel::updateJ()
    {
        volVectorField drivingForce =
        (
            fvc::div(SigmaR()) + divSigmaL1() + divSigmaL2() + Fp()
        );

        dimensionedScalar dragCoefficient = ((2.0*sqr(a_))/(9.0*nu0_));

        J_ = dragCoefficient*drivingForce*hindranceModelPtr_->fH();
        J_.correctBoundaryConditions();
    }


    void Foam::suspensionModel::updateDivSigmaL2f()
    {
        divSigmaL2f_ =
        (
            -fvc::interpolate(nun())*(fvc::snGrad(gammaDot())*Qf()
            + fvc::snGrad(gammaNL())*mesh_.magSf())
            -fvc::interpolate(nun()*gammaDot())*divQf()
        );
    }

} // End namespace Foam

// ************************************************************************* //

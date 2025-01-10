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

Class
    Foam::suspensionModel

Description
    Base class for creation of suspension models.

SourceFiles
    suspensionModel.C

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

#ifndef suspensionModel_H
#define suspensionModel_H

#include "volFields.H"
#include "surfaceInterpolate.H"
#include "IOdictionary.H"
#include "fvCFD.H"
#include "nonLocalShearRateModel.H"
#include "anisotropicStressTensorModel.H"
#include "shearViscosityModel.H"
#include "normalViscosityModel.H"
#include "hindranceModel.H"
#include "solidsPressureModel.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class suspensionModel Declaration
\*---------------------------------------------------------------------------*/

class suspensionModel : public IOdictionary
{

    private:

    //- Private data
    const fvMesh &mesh_;
    dictionary suspensionDict_;
    const volScalarField& p_;
    const volVectorField& U_;
    dimensionedScalar a_;
    dimensionedScalar nu0_;
    dimensionedScalar cMax_;
    volVectorField J_;
    volTensorField Sigma_;
    volScalarField gamma_;
    volScalarField gammaNL_;
    const volScalarField& c_;
    const surfaceScalarField& phi_;
    const dimensionedTensor I_;
    volVectorField divSigmaL1_;
    surfaceScalarField divSigmaL1f_;
    surfaceScalarField divSigmaL2f_;
    autoPtr<nonLocalShearRateModel> nonLocalShearRateModelPtr_;
    autoPtr<anisotropicStressTensorModel> anisotropicStressTensorModelPtr_;
    autoPtr<shearViscosityModel> shearViscosityModelPtr_;
    autoPtr<normalViscosityModel> normalViscosityModelPtr_;
    autoPtr<hindranceModel> hindranceModelPtr_;
    autoPtr<solidsPressureModel> solidsPressureModelPtr_;


    //- Disallow copy construct
    suspensionModel(const suspensionModel&);


    //- Disallow default bitwise assignment
    void operator=(const suspensionModel&);

    public:

    //- Construct from components
    suspensionModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const volScalarField& p,
        const volVectorField& U,
        const volScalarField& alpha,
        const surfaceScalarField& phi
    );


    //- Destructor

    virtual ~suspensionModel()
    {}


    //- Member functions
    void update();
    void updateGamma();
    void updateDivSigmaL1f();
    void updateJ();
    void updateDivSigmaL2f();
    const volScalarField& p() const
    {
        return p_;
    }
    const volVectorField& U() const
    {
        return U_;
    }
    const volScalarField& c() const
    {
        return c_;
    }
    const surfaceScalarField& phi() const
    {
        return phi_;
    }
    const dimensionedTensor I() const
    {
        return I_;
    }
    tmp<volTensorField> Q() const
    {
        return anisotropicStressTensorModelPtr_->Q();
    }
    tmp<surfaceScalarField> Qf() const
    {
        return fvc::interpolate(Q()) & (mesh_.Sf()/mesh_.magSf()) & mesh_.Sf();
    }
    tmp<surfaceScalarField> divQf() const
    {
        return (fvc::interpolate(fvc::div(Q())) & mesh_.Sf());
    }
    dimensionedScalar nu0() const
    {
        return nu0_;
    }
    tmp<volScalarField> nus() const
    {
        return shearViscosityModelPtr_->nus();
    }
    tmp<volScalarField> nup() const
    {
        return shearViscosityModelPtr_->nup();
    }
    tmp<volScalarField> nun() const
    {
        return normalViscosityModelPtr_->nun();
    }
    tmp<volScalarField> dNunDc() const
    {
        return normalViscosityModelPtr_->dNunDc();
    }
    tmp<volTensorField> GrU() const
    {
        return fvc::grad(U_);
    }
    tmp<volTensorField> E() const
    {
        return ((symm(GrU())) & I());
    }
    tmp<volScalarField> gammaDot() const
    {
        return sqrt((2.0*(E() && E())));
    }
    tmp<volScalarField> gammaNL() const
    {
        return nonLocalShearRateModelPtr_->gammaNL();
    }
    tmp<volTensorField> SigmaR() const
    {
        return 2.0*nup()*E();
    }
    tmp<volTensorField> SigmaF() const
    {
        return (-p_*I() + 2.0*nu0_*E());
    }
    tmp<volTensorField> SigmaL() const
    {
        return -nun()*(gammaDot()*Q() + gammaNL()*I());
    }
    tmp<volTensorField> SigmaP() const
    {
        return (SigmaL() + SigmaR());
    }
    tmp<volTensorField> Sigma()
    {
        Sigma_ = SigmaP() + SigmaF();
        return Sigma_;
    }
    tmp<volVectorField> divSigmaL1() const
    {
        return divSigmaL1_;
    }
    tmp<surfaceScalarField> divSigmaL1f() const
    {
        return divSigmaL1f_;
    }
    tmp<volVectorField> divSigmaL2() const
    {
        return -nun()*(fvc::div(gammaDot()*Q()) + fvc::div(gammaNL()*I()));
    }
    tmp<surfaceScalarField> divSigmaL2f() const
    {
        return divSigmaL2f_;
    }
    tmp<volVectorField> Fp() const
    {
        return solidsPressureModelPtr_->Fp();
    }
    tmp<surfaceScalarField> Fpf() const
    {
        return solidsPressureModelPtr_->Fpf();
    }
    tmp<volVectorField> J()
    {
        return J_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

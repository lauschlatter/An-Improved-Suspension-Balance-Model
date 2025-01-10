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

InClass
    Foam::improvedFrameInvariantAnisotropicStressTensor

\*---------------------------------------------------------------------------*/

#include "improvedFrameInvariantAnisotropicStressTensor.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyFvsPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
    defineTypeNameAndDebug(improvedFrameInvariantAnisotropicStressTensor, 0);

    addToRunTimeSelectionTable
    (
        anisotropicStressTensorModel,
        improvedFrameInvariantAnisotropicStressTensor,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicStressTensorModels::
improvedFrameInvariantAnisotropicStressTensor::
improvedFrameInvariantAnisotropicStressTensor
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    anisotropicStressTensorModel(mesh, U, c, dict),
    improvedFrameInvariantAnisotropicStressTensorDict_
    (
        dict.subDict(typeName+"Coeffs")
    ),
    lambda1_
    (
        improvedFrameInvariantAnisotropicStressTensorDict_.lookup("lambda1")
    ),
    lambda2_
    (
        improvedFrameInvariantAnisotropicStressTensorDict_.lookup("lambda2")
    ),
    lambda3_
    (
        improvedFrameInvariantAnisotropicStressTensorDict_.lookup("lambda3")
    ),
    cut_(improvedFrameInvariantAnisotropicStressTensorDict_.lookup("cut")),
    decay_(improvedFrameInvariantAnisotropicStressTensorDict_.lookup("decay")),
    rTol_(improvedFrameInvariantAnisotropicStressTensorDict_.lookup("rTol")),
    gUCrit_("gUCrit", dimensionSet(0,0,0,0,0,0,0), 0.0),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("epsilon", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "calculated"
    ),
    fx_
    (
        IOobject
        (
            "fx_",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("fx_", dimensionSet(0,0,0,0,0,0,0), scalar(0)),
        "calculated"
    ),
    y_
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("y", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "calculated"
    ),
    dirU_
    (
        IOobject
        (
            "dirU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE//IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("dirU", dimensionSet(0,0,0,0,0,0,0), vector(0,0,0)),
        "zeroGradient"
    ),
    dirGradU_
    (
        IOobject
        (
            "dirGradU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE//IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "dirGradU", dimensionSet(0,0,0,0,0,0,0), vector(0,0,0)
        ),
        "zeroGradient"
    ),
    dirCurlU_
    (
        IOobject
        (
            "dirCurlU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE//IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "dirCurlU", dimensionSet(0,0,0,0,0,0,0), vector(0,0,0)
        ),
        "zeroGradient"
    ),
    curlU_
    (
        IOobject
        (
            "curlU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("curlU", dimensionSet(0,0,-1,0,0,0,0), vector(0,0,0)),
        "zeroGradient"
    ),
    gradU_
    (
        IOobject
        (
            "gradU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "gradU", dimensionSet(0,0,-1,0,0,0,0), tensor(0,0,0,0,0,0,0,0,0)
        ),
        "zeroGradient"
    ),
    magGradU_
    (
        IOobject
        (
            "magGradU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("magGradU", dimensionSet(0,0,-1,0,0,0,0), scalar(0)),
        "zeroGradient"
    ),
    magCurlU_
    (
        IOobject
        (
            "magCurlU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("magCurlU", dimensionSet(0,0,-1,0,0,0,0), scalar(0)),
        "zeroGradient"
    ),
    magU_
    (
        IOobject
        (
            "magU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("magU", dimensionSet(0,1,-1,0,0,0,0), scalar(0)),
        "zeroGradient"
    ),
    l1_
    (
        IOobject
        (
            "l1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("l1", dimensionSet(0,0,0,0,0,0,0), lambda1_.value())
    ),
    l2_
    (
        IOobject
        (
            "l2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("l2", dimensionSet(0,0,0,0,0,0,0), lambda2_.value())
    ),
    l3_
    (
        IOobject
        (
            "l3",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("l3", dimensionSet(0,0,0,0,0,0,0), lambda3_.value())
    ),
    is2D_(false),
    e32D_
    (
        improvedFrameInvariantAnisotropicStressTensorDict_.lookupOrDefault
        (
            "e3",
            vector(0,0,1)
        )
    ),
    useGradients_
    (
        improvedFrameInvariantAnisotropicStressTensorDict_.lookupOrDefault
        (
            "useGradients",
            false
        )
    )
{
    wordList UTypes
    (
        U.boundaryField().size()
    );

    Info<< "Checking boundary field" << endl;
    forAll(U.boundaryField(), i)
    {
        UTypes[i] = U.boundaryField()[i].type() ;
        if (UTypes[i] == "empty")
        {
            is2D_ = true;
            Info<< "Found an empty boundary condition in patch." << endl;
            Info<< "This case will be treated as 2D." << endl;
            Info<< "e3 direction will be set as " << e32D_ << endl;
            dirCurlU_ = e32D_;
            dirCurlU_.correctBoundaryConditions();
        }
    }

    dimensionedScalar yHat = (1/decay_)*log(rTol_/(1 - rTol_));
    gUCrit_ = cut_.value() - tanh(yHat);
    //Info<< "\tcalculated gUCrit:" << gUCrit_.value() << endl;
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::anisotropicStressTensorModels::
improvedFrameInvariantAnisotropicStressTensor::
updateAnisotropicStressTensor()
{
    dimensionedVector i ("i", dimensionSet(0,0,-1,0,0,0,0), vector(1,0,0));
    dimensionedVector j ("j", dimensionSet(0,0,-1,0,0,0,0), vector(0,1,0));
    dimensionedVector k ("k", dimensionSet(0,0,-1,0,0,0,0), vector(0,0,1));

    curlU_ = fvc::curl(U_);
    curlU_.correctBoundaryConditions();
    magCurlU_ = mag(curlU_);
    magCurlU_.correctBoundaryConditions();
    magU_ = mag(U_);
    magU_.correctBoundaryConditions();
    gradU_ = fvc::grad(U_);
    gradU_.correctBoundaryConditions();
    magGradU_ = sqrt(gradU_ && gradU_);
    magGradU_.correctBoundaryConditions();

    dimensionedScalar maxMagGradU = max(magGradU_);

    dimensionedScalar admPar =
        (
            max(magU_)/max
            (
                max(magCurlU_),
                dimensionedScalar("curlV", dimVelocity/dimLength, 1e-18)
            )
        );

    // e_1
    dirU_ = U_/(max(mag(U_), dimensionedScalar("V", dimVelocity, 1e-12)));
    dirU_.correctBoundaryConditions();

    if (is2D_)
    {}
    else
    {
        // e_3
        dirCurlU_ =
            (
                fvc::curl(U_)/max
                    (
                        mag(fvc::curl(U_)),
                        dimensionedScalar("curlV", dimVelocity/dimLength, 1e-12)
                    )
            );

        dirCurlU_.correctBoundaryConditions();

        volVectorField ax = dirU_ ^ i;
        volVectorField ay = dirU_ ^ j;
        volVectorField az = dirU_ ^ k;
        volVectorField aa = ax;

        if (useGradients_)
        {
            epsilon_ = magGradU_/max
                (
                    maxMagGradU,
                    dimensionedScalar("magGradU", dimVelocity/dimLength, 1e-12)
                );
        }
        else
        {
            epsilon_ = magCurlU_*admPar.value()/max
                (
                    magU_,
                    dimensionedScalar("magGradU", dimVelocity/dimLength, 1e-12)
                );
        }

        epsilon_.correctBoundaryConditions();

        forAll(magCurlU_, cellI)
        {
            if (epsilon_[cellI] < cut_.value())
            {
                if (mag(ax[cellI]) > mag(ay[cellI]))
                {
                    if (mag(ax[cellI]) > mag(az[cellI]))
                    {
                        aa[cellI] = ax[cellI];
                    }
                    else
                    {
                        aa[cellI] = az[cellI];
                    }
                }
                else
                {
                    if (mag(ay[cellI]) > mag(az[cellI]))
                    {
                        aa[cellI] = ay[cellI];
                    }
                    else
                    {
                        aa[cellI] = az[cellI];
                    }
                }
                dirCurlU_[cellI] = (aa[cellI])/max
                        (
                            mag(aa[cellI]),small
                        );
            }
        }
    y_ = atanh(epsilon_ - gUCrit_);
    fx_ = pow(max(1. + exp(decay_*y_), 1e-12), -1.);
    l3_ = fx_*(lambda2_ - lambda3_) + lambda3_;
    Info<< "lambda3 Max = " << max(l3_).value()
        << "; lambda3 Min = " << min(l3_).value() << endl;
    }

    // e_2
    dirGradU_ = dirU_ ^ dirCurlU_/max
        (
            mag(dirU_ ^ dirCurlU_),
            dimensionedScalar("curlV", dimless, 1e-12)
        );
    dirGradU_.correctBoundaryConditions();
    

    volTensorField Qrot =
        (
            l1_*(dirU_*dirU_)
            + l2_*(dirGradU_*dirGradU_)
            + l3_*(dirCurlU_*dirCurlU_)
        );

    Q_ = Qrot;
    Q_.correctBoundaryConditions();

}

// -------------------------------------------------------------------------//

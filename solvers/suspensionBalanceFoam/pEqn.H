p.boundaryFieldRef().updateCoeffs();
volScalarField rUA = 1.0/UEqn.A();
surfaceScalarField rUAf = fvc::interpolate(rUA);
U = rUA*UEqn.H();
tUEqn.clear();
phi =
(
    (fvc::interpolate(U) & mesh.Sf())
    + mixture.divSigmaL1f()*rUAf + mixture.Fpf()*rUAf
);

adjustPhi(phi, U, p);

// Store pressure for under-relaxation
p.storePrevIter();

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUA, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();

    if (piso.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
U += (mixture.divSigmaL1() + mixture.Fp() - fvc::grad(p))*rUA;
U.correctBoundaryConditions();

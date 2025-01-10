mixture.updateJ();

// Construct and solve the particle continuity equation
fvScalarMatrix cEqn
(
    fvm::ddt(c) + fvm::div(phi, c) == - fvc::div(mixture.J())
);

cEqn.solve();
c.correctBoundaryConditions();

Info<< "min(c) = " << min(c).value()
    << "\t ave(c) = " << c.weightedAverage(mesh.V()).value()
    << "\t max(c) = " << max(c).value()
    << endl;

if (fixCAve)
{
    Info<< "Fixing CAve values based on volume average" << endl;

    c = c*cAve/c.weightedAverage(mesh.V()).value();

    Info<< "min(c) = " << min(c).value()
        << "\t ave(c) = " << c.weightedAverage(mesh.V()).value()
        << "\t max(c) = " << max(c).value()
        << endl;
}

if (fixCflux)
{
    Info<< "Fixing CAve values based on inlet flux" << endl;

    scalar ain(0);
    scalar volflux(0);
    scalar influx(0);

    const scalarField& aIn = mesh.magSf().boundaryField()[inletPatchID];
    const scalarField& cIn = c.boundaryField()[inletPatchID];
    const scalarField& phiIn = phi.boundaryField()[inletPatchID];

    forAll(cIn, faceI)
    {
        volflux += phiIn[faceI];
        influx += phiIn[faceI]*cIn[faceI];
        ain += aIn[faceI];
    }

    reduce(volflux, sumOp<scalar>());
    reduce(influx, sumOp<scalar>());
    reduce(ain, sumOp<scalar>());

    c = c*cAve*volflux/(influx);

    Info<< "min(c) = " << min(c).value()
        << "\t ave(c) = " << c.weightedAverage(mesh.V()).value()
        << "\t max(c) = " << max(c).value()
        << endl;

    Info<< "Average Uin = " << -volflux/ain << endl;
    Info<< "Set point flux = " << cAve*volflux << endl;
    Info<< "Calculated flux = " << influx << endl;
}

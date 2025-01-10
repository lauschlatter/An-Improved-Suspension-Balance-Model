// Update flux before momentum predictor
mixture.updateDivSigmaL1f();

// Momentum predictor
tmp<fvVectorMatrix> tUEqn
(
    - fvm::laplacian(mixture.nus(), U)
    - fvc::div(mixture.nus()*dev(fvc::grad(U)().T()))
    - mixture.divSigmaL2()
    ==
    fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

fvOptions.constrain(UEqn);

UEqn.relax();

if (piso.momentumPredictor())
{
    solve(UEqn == mixture.divSigmaL1() + mixture.Fp() - fvc::grad(p));
    fvOptions.correct(U);
}

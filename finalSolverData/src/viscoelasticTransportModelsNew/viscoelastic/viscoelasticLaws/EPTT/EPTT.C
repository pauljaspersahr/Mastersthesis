/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "EPTT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EPTT, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, EPTT, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EPTT::EPTT
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const areaVectorField& Us, // added
    const edgeScalarField& phis, // added
    const dictionary& dict
)
:   
    viscoelasticLaw(name, U, phi, Us, phis), // changed
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    tauSurface_ // added
    (
        IOobject
        (
            "tauSurface" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Us.mesh()
    ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    epsilon_(dict.lookup("epsilon")),
    lambda_(dict.lookup("lambda")),
    zeta_(dict.lookup("zeta"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// added
void Foam::EPTT::tauSurfaceCorrect()
{

    const tmp<areaTensorField> tL = fac::grad(Us());
    const areaTensorField& L = tL();

    // Convected derivate term
    areaTensorField C = tauSurface_ & L;

    // Twice the rate of deformation tensor
    areaSymmTensorField twoD = twoSymm(L);

     // Stress transport equation
    Info<< "creating faSymmTensorMatrix"<< endl;

    faSymmTensorMatrix tauSurfaceEqn
    (
        lambda_*fam::ddt(tauSurface_)
      + lambda_*fam::div(phis(), tauSurface_)
     ==
        etaP_*twoD
      + lambda_*twoSymm(C)
      - lambda_*zeta_*symm(tauSurface_ & twoD)
      - fam::Sp
        (
            Foam::exp(epsilon_*lambda_/etaP_*tr(tauSurface_)),
            tauSurface_
        )

    );

    tauSurfaceEqn.relax();
    tauSurfaceEqn.solve(Us().mesh().solutionDict().subDict(tauSurface_.name()));


}



Foam::tmp<Foam::fvVectorMatrix> Foam::EPTT::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_, "div(tau)")
      - fvc::laplacian(etaPEff, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_), U, "laplacian(etaPEff+etaS,U)")
    );
}


void Foam::EPTT::tauCorrect()
{
    // Velocity gradient tensor
    const tmp<volTensorField> tL = fvc::grad(U());
    const volTensorField& L = tL();

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

     // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        lambda_*fvm::ddt(tau_)
      + lambda_*fvm::div(phi(), tau_)
     ==
        etaP_*twoD
      + lambda_*twoSymm(C)
      - lambda_*zeta_*symm(tau_ & twoD)
      - fvm::Sp
        (
            Foam::exp(epsilon_*lambda_/etaP_*tr(tau_)),
            tau_
        )
    );

    tauEqn.relax();
    tauEqn.solve();
}


// ************************************************************************* //

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

#include "viscoelasticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(viscoelasticModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

viscoelasticModel::viscoelasticModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const areaVectorField& Us,      // added
    const edgeScalarField& phis     // added
)
:
    IOdictionary
    (
        IOobject
        (
            "viscoelasticProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    lawPtr_(viscoelasticLaw::New(word::null, U, phi, Us, phis, subDict("rheology")))    // edited
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// added
tmp<areaSymmTensorField> viscoelasticModel::tauSurface() const
{
    return lawPtr_->tauSurface();
}

// added
void viscoelasticModel::tauSurfaceCorrect()
{
    lawPtr_->tauSurfaceCorrect();
}

tmp<volSymmTensorField> viscoelasticModel::tau() const
{
    return lawPtr_->tau();
}


tmp<fvVectorMatrix> viscoelasticModel::divTau(volVectorField& U) const
{
    return lawPtr_->divTau(U);
}

// edited
void viscoelasticModel::tauCorrect()
{
    lawPtr_->tauCorrect();
}


bool viscoelasticModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

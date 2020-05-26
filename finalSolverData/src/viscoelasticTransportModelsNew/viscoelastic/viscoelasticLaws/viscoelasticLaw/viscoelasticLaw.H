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

Class
    viscoelasticLaw

Description
    Abstract base class for viscoelastic flow models.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    viscoelasticLaw.C
    newViscoelasticLaw.C

\*---------------------------------------------------------------------------*/

#ifndef viscoelasticLaw_H
#define viscoelasticLaw_H

#include "volFields.H"
#include "surfaceFields.H"
#include "areaFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// added: finite area method
#include "faCFD.H"                                                                 
#include "fam.H"
#include "fac.H"
#include "faMatrices.H"
#include "faMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class viscoelasticLaw Declaration
\*---------------------------------------------------------------------------*/

class viscoelasticLaw
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to face flux field
        const surfaceScalarField& phi_;

        //- added: Reference to surface velocity field
        const areaVectorField& Us_;

        //- added: Reference to edge flux field
        const edgeScalarField& phis_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscoelasticLaw(const viscoelasticLaw&);

        //- Disallow default bitwise assignment
        void operator=(const viscoelasticLaw&);


public:

    //- Runtime type information
    TypeName("viscoelasticLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            viscoelasticLaw,
            dictionary,
            (
                const word& name,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const areaVectorField& Us,      // added
                const edgeScalarField& phis,    // added
                const dictionary& dict
            ),
            (name, U, phi, Us, phis, dict)      // edited
        );


    // Constructors

        //- Construct from components
        viscoelasticLaw
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const areaVectorField& Us,          // added
            const edgeScalarField& phis         // added
        );


    // Selectors

        //- Return a reference to the selected viscoelastic law
        static autoPtr<viscoelasticLaw> New
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const areaVectorField& Us,          // added
            const edgeScalarField& phis,        // added
            const dictionary& dict
        );


    // Destructor

        virtual ~viscoelasticLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the face flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- added: Return the surface velocity field
        const areaVectorField& Us() const       
        {
            return Us_;
        }

        //- added: Return the edge flux field
        const edgeScalarField& phis() const
        {
            return phis_;                       
        }

        //- added: Return viscoelastic surface stress tensor
        virtual tmp<areaSymmTensorField> tauSurface() const = 0;

        //- added: Correct viscoelastic surface stress
        virtual void tauSurfaceCorrect() = 0;

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const = 0;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const = 0;

        //- edited: Correct the viscoelastic stress
        virtual void tauCorrect() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

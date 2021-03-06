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
    Foam::fvPatch

Description
    Foam::fvPatch

SourceFiles
    fvPatch.C
    newFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatch_H
#define fvPatch_H

#include "polyPatch.H"
#include "labelList.H"
#include "SubList.H"
#include "typeInfo.H"
#include "tmp.H"
#include "primitiveFields.H"
#include "SubField.H"
#include "fvPatchFieldsFwd.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvBoundaryMesh;
class surfaceInterpolation;

/*---------------------------------------------------------------------------*\
                           Class fvPatch Declaration
\*---------------------------------------------------------------------------*/

class fvPatch
{
    // Private data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;

        //- Reference to boundary mesh
        const fvBoundaryMesh& boundaryMesh_;


    // Private Member Functions

        //- Disallow construct as copy
        fvPatch(const fvPatch&);

        //- Disallow assignment
        void operator=(const fvPatch&);


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const;

        //- Make patch face non-orthogonality correction vectors
        virtual void makeCorrVecs(vectorField&) const;

        //- Initialise the patches for moving points
        virtual void initMovePoints();

        //- Correct patches after moving points
        virtual void movePoints();


public:

    typedef fvBoundaryMesh BoundaryMesh;

    friend class fvBoundaryMesh;
    friend class surfaceInterpolation;

    //- Runtime type information
    TypeName(polyPatch::typeName_());


    // Declare run-time constructor selection tables

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            fvPatch,
            polyPatch,
            (const polyPatch& patch, const fvBoundaryMesh& bm),
            (patch, bm)
        );
#endif


    // Constructors

        //- Construct from polyPatch
        fvPatch(const polyPatch&, const fvBoundaryMesh&);


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<fvPatch> New(const polyPatch&, const fvBoundaryMesh&);


    // Destructor

        virtual ~fvPatch();


    // Member Functions

    // Access

            //- Return the polyPatch
            const polyPatch& patch() const
            {
                return polyPatch_;
            }

            //- Return name
            const word& name() const
            {
                return polyPatch_.name();
            }

            //- Return size
            virtual label size() const
            {
                return polyPatch_.size();
            }

            //- Return true if this patch is coupled
            bool coupled() const
            {
                return polyPatch_.coupled();
            }

            //- Return true if this patch is wall
            bool isWall() const
            {
                return polyPatch_.isWall();
            }

            //- Return true if the given type is a constraint type
            static bool constraintType(const word& pt);

            //- Return a list of all the constraint patch types
            static wordList constraintTypes();

            //- Return the index of this patch in the fvBoundaryMesh
            label index() const
            {
                return polyPatch_.index();
            }

            //- Return boundaryMesh reference
            const fvBoundaryMesh& boundaryMesh() const
            {
                return boundaryMesh_;
            }

            //- Slice list to patch
            template<class T>
            const typename List<T>::subList patchSlice(const List<T>& l) const
            {
                return typename List<T>::subList
                (
                    l, size(),
                    polyPatch_.start()
                );
            }

            template<class T>
            const typename Field<T>::subField patchSlice
            (
                const Field<T>& l
            ) const
            {
                return typename Field<T>::subField
                (
                    l, size(),
                    polyPatch_.start()
                );
            }

            //- Return faceCells
            virtual const unallocLabelList& faceCells() const;


        // Access functions for geometrical data

            //- Return face centres
            const vectorField& Cf() const;

            //- Return neighbour cell centres
            tmp<vectorField> Cn() const;

            //- Return face area vectors
            const vectorField& Sf() const;

            //- Return face area magnitudes
            const scalarField& magSf() const;

            //- Return face normals
            tmp<vectorField> nf() const;

            //- Return cell-centre to face-centre vector
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre vector is returned
            virtual tmp<vectorField> delta() const;


        // Access functions for demand driven data

            //- Return patch weighting factors
            const scalarField& weights() const;

            //- Return the face - cell distance coeffient
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre distance coeffient is returned
            const scalarField& deltaCoeffs() const;


        // Evaluation functions

            //- Return given internal field next to patch as patch field
            template<class Type>
            tmp<Field<Type> > patchInternalField(const UList<Type>&) const;

            //- Return the corresponding patchField of the named field
            template<class GeometricField, class Type>
            const typename GeometricField::PatchFieldType& patchField
            (
                const GeometricField&
            ) const;

            //- Lookup and return the patchField of the named field from the
            //  local objectRegistry.
            //  N.B.  The dummy pointer arguments are used if this function is
            //  instantiated within a templated function to avoid a bug in gcc.
            //  See inletOutletFvPatchField.C and outletInletFvPatchField.C
            template<class GeometricField, class Type>
            const typename GeometricField::PatchFieldType& lookupPatchField
            (
                const word& name,
                const GeometricField* = NULL,
                const Type* = NULL
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

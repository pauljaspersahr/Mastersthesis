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

Description

\*---------------------------------------------------------------------------*/

#include "freeSurface.H"

#include "volFields.H"
#include "transformField.H"

#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(modelPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurface::freeSurface
(
    dynamicFvMesh& m,
    const volScalarField& rho,
    volVectorField& Ub,
    volScalarField& Pb,
    const surfaceScalarField& sfPhi                                                                   
)
:
    IOdictionary
    (
        IOobject
        (
            "freeSurfaceProperties",
            Ub.mesh().time().constant(),
            Ub.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(m),
    rho_(rho),
    U_(Ub),
    p_(Pb),
    phi_(sfPhi),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_
    (
        this->lookup("twoFluids")
    ),
    normalMotionDir_
    (
        this->lookup("normalMotionDir")
    ),
    motionDir_(0, 0, 0),
    cleanInterface_
    (
        this->lookup("cleanInterface")
    ),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_
    (
        "muFluidA",dimensionSet(1,-1,-1,0,0,0,0), scalar(0)       // used for FS Velocity update as total viscosity
    ),
    muFluidB_
    (
        this->lookup("muFluidB")                            // inside, mu_s solvent contribution
    ),
    rhoFluidA_
    (
        this->lookup("rhoFluidA")
    ),
    rhoFluidB_
    (
        this->lookup("rhoFluidB")
    ),
    g_(this->lookup("g")),
    cleanInterfaceSurfTension_
    (
        this->lookup("surfaceTension")
    ),
    fixedFreeSurfacePatches_
    (
        this->lookup("fixedFreeSurfacePatches")
    ),
    pointNormalsCorrectionPatches_
    (
        this->lookup("pointNormalsCorrectionPatches")
    ),
    nFreeSurfCorr_
    (
        readInt(this->lookup("nFreeSurfaceCorrectors"))
    ),
    smoothing_(false),
    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    totalDisplacementPtr_(NULL),
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL),
    modelPtr_(NULL)
//     lawPtrPtr_(NULL)
    //lawPtr_(viscoelasticLaw::New(word::null, Ub, sfPhi, UsPtr_, phisPtr_, subDict("rheology")))                    // PJS, added
{


	Info << "Calculating muFluidA" << endl;
    
    IOdictionary viscoelasticProperties
    (
        IOobject
        (
            "viscoelasticProperties",
            Ub.mesh().time().constant(),
            Ub.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //Sum up solvent related viscosities
    PtrList<entry> modelEntries_(viscoelasticProperties.subDict("rheology").lookup("models"));                         // added PJS, 19.08.19

	Info << "Model entry list created" << endl;

    forAll(modelEntries_, modelI)
    {
        muFluidA_ += dimensionedScalar(modelEntries_[modelI].dict().lookup("etaS")); 
    }
    Info << "muFluidA cualculated as:" << muFluidA_.value() << endl;





    //Read motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = vector(this->lookup("motionDir"));
        motionDir_ /= mag(motionDir_) + SMALL;
    }

    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "freeSurface::freeSurface(...)"
            )   << "Patch name for point normals correction does not exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }

    // Clear geometry
    aMesh().movePoints();


    // Detect the free surface patch
    forAll (mesh().boundary(), patchI)
    {
        if(mesh().boundary()[patchI].name() == "freeSurface")
        {
            aPatchID_ = patchI;

            Info<< "Found free surface patch. ID: " << aPatchID_
                << endl;
        }
    }

    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Free surface patch not defined.  Please make sure that "
                << " the free surface patches is named as freeSurface"
                << abort(FatalError);
    }


    // Detect the free surface shadow patch
    if (twoFluids())
    {
        forAll (mesh().boundary(), patchI)
        {
            if(mesh().boundary()[patchI].name() == "freeSurfaceShadow")
            {
                bPatchID_ = patchI;

                Info<< "Found free surface shadow patch. ID: "
                    << bPatchID_ << endl;
            }
        }

        if(bPatchID() == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Free surface shadow patch not defined. "
                    << "Please make sure that the free surface shadow patch "
                    << "is named as freeSurfaceShadow."
                    << abort(FatalError);
        }
    }


    // Mark free surface boundary points
    // which belonge to processor patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }
        }
    }


    // Mark fixed free surface boundary points
    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }


    // Mark free-surface boundary point
    // at the axis of 2-D axisymmetic cases
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            if(wedgePatch.axisPoint() > -1)
            {
                motionPointsMask()[wedgePatch.axisPoint()] = 0;

                Info << "Axis point: "
                    << wedgePatch.axisPoint()
                    << "vector: "
                    << aMesh().points()[wedgePatch.axisPoint()] << endl;
            }
        }
    }


    // Read free-surface points total displacement if present
    readTotalDisplacement();


    // Read control points positions if present
    controlPoints();


    // Check if smoothing switch is set
    if (this->found("smoothing"))
    {
        smoothing_ = Switch(this->lookup("smoothing"));
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

freeSurface::~freeSurface()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //






volSymmTensorField freeSurface::calculateTauTotalField(volVectorField& U_, volScalarField& fluidIndicator_) const
{

    
	// Velocity gradient tensor
    const tmp<volTensorField> tL = fvc::grad(U_);
    const volTensorField& L = tL();


    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

	volSymmTensorField tau_ = (fluidIndicator_ * (muFluidA_ - muFluidB_) + muFluidB_ ) * twoD;

	// tau_ += fluidIndicator_ * tau()();

	return tau_;
}


tmp<areaSymmTensorField> freeSurface::tauSurface() const
{
    return modelPtr().tauSurface();
}

void freeSurface::tauSurfaceCorrect()
{
    modelPtr().tauSurfaceCorrect();
}


tmp<volSymmTensorField> freeSurface::tau()                                       // added
{
    return modelPtr().tau();
}


tmp<fvVectorMatrix> freeSurface::divTau(volVectorField& U_)                       // added
{
    return modelPtr().divTau(U_);
}


void freeSurface::tauCorrect()                                                           // added
{
    modelPtr().tauCorrect();
}


bool freeSurface::read()                                                              // added
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



void freeSurface::updateDisplacementDirections()
{
    if(normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correcte point displacement direction
        // at the "centerline" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll(aMesh().boundary(), patchI)                                          // loops over all cells from aMesh.boundary till patchI
        {
            if(aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

                vector axis = wedgePatch.axis();

                label centerLinePatchID =
                    aMesh().boundary().findPatchID("centerline");

                if(centerLinePatchID != -1)
                {
                    const labelList& pointLabels =
                        aMesh().boundary()[centerLinePatchID].pointLabels();

                    forAll(pointLabels, pointI)
                    {
                        vector dir =
                            pointsDisplacementDir()[pointLabels[pointI]];

                        dir = (dir&axis)*axis;
                        dir /= mag(dir);

                        pointsDisplacementDir()[pointLabels[pointI]] = dir;
                    }
                }
                else
                {
                    Info << "Warning: centerline polyPatch does not exist. "
                        << "Free surface points displacement directions "
                        << "will not be corrected at the axis (centerline)"
                        << endl;
                }

                break;
            }
        }

        // Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() =
            facesDisplacementDir()
           *(facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }
}


bool freeSurface::predictPoints()
{
    // Smooth interface

    if (smoothing_)
    {
        controlPoints() = aMesh().areaCentres().internalField();
        movePoints(scalarField(controlPoints().size(), 0));
        movePoints(-fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()]);
    }

    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool freeSurface::correctPoints()
{
    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
    }

    return true;
}


bool freeSurface::movePoints(const scalarField& interfacePhi)                           // interfacephi == velocity on the free surface?
{
    pointField newMeshPoints = mesh().points();                                         // write pointField from mesh points

    scalarField sweptVolCorr =                                                          // sweptVolCorr == ?
        interfacePhi
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];                          // meshPhi: Calculates the mesh motion flux and converts fluxes from absolute to relative and back

    word ddtScheme
    (
        mesh().schemesDict().ddtScheme
        (
            "ddt(" + rho().name() + ',' + U().name()+')'
        )
    );

    if
    (
        ddtScheme
     == fv::CrankNicolsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
    }
    else
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    const scalarField& Sf = aMesh().S();                                    // Sf == Surface?
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();      // normals

    scalarField deltaH =                                                    // deltaH ? 
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    pointField displacement = pointDisplacement(deltaH);


    // Move only free-surface points

    const labelList& meshPointsA =
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB =
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );

        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI];
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    // faMesh motion is done automatically, using meshObject
    // HJ, 8/Aug/2011
//     aMesh().movePoints(mesh().points());


    // Move correctedFvPatchField fvSubMeshes                                           // 

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }


        //Update tau boundary conditions, added 27.08. PJS
        forAll(tau()().boundaryField(), patchI)
        {
            if
            (
                (
                    tau()().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<symmTensor>::typeName
                )
                ||
                (
                    tau()().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<symmTensor>::typeName
                )
                ||
                (
                    tau()().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<symmTensor>::typeName
                )
            )
            {
                correctedFvPatchField<symmTensor>& aTau =
                    refCast<correctedFvPatchField<symmTensor> >
                    (
                        tau()().boundaryField()[patchI]
                    );

                aTau.movePatchSubMesh();
            }
        }


    return true;
}


bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
        pointField newPoints = mesh().points();

        const labelList& meshPointsA =
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        forAll (totalDisplacement(), pointI)
        {
            newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI];               // a -=  b     ==      a = a - b
        }


        // Check mesh motion solver type
        bool feMotionSolver =
            mesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );
        bool fvMotionSolver =
            mesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUaPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            tetPolyPatchInterpolation tppiAPatch
            (
                refCast<const faceTetPolyPatch>
                (
                    motionUaPatch.patch()
                )
            );

            motionUaPatch ==
                tppiAPatch.pointToPointInterpolate
                (
                    totalDisplacement()/DB().deltaT().value()
                );

            if(twoFluids_)
            {
                const labelList& meshPointsB =
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );

                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -=
                        totDisplacementB[pointI];
                }

                fixedValueTetPolyPatchVectorField& motionUbPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                tetPolyPatchInterpolation tppiBPatch
                (
                    refCast<const faceTetPolyPatch>(motionUbPatch.patch())
                );

                motionUbPatch ==
                    tppiBPatch.pointToPointInterpolate
                    (
                        totDisplacementB/DB().deltaT().value()
                    );
            }
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    mesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUaPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[aPatchID()]
                );

            motionUaPatch ==
                totalDisplacement()/DB().deltaT().value();

            if(twoFluids_)
            {
                const labelList& meshPointsB =
                    mesh().boundaryMesh()[bPatchID()].meshPoints();

                pointField totDisplacementB =
                    interpolatorAB().pointInterpolate
                    (
                        totalDisplacement()
                    );

                forAll (totDisplacementB, pointI)
                {
                    newPoints[meshPointsB[pointI]] -=
                        totDisplacementB[pointI];
                }

                fixedValuePointPatchVectorField& motionUbPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[bPatchID()]
                    );

                motionUbPatch ==
                    totDisplacementB/DB().deltaT().value();
            }
        }

        twoDPointCorrector twoDPointCorr(mesh());

        twoDPointCorr.correctPoints(newPoints);

        mesh().movePoints(newPoints);

        deleteDemandDrivenData(totalDisplacementPtr_);

        mesh().update();

        // faMesh motion is done automatically, using meshObject
        // HJ, 8/Aug/2011
//         aMesh().movePoints(mesh().points());

        // Move correctedFvPatchField fvSubMeshes                                       /

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }
    

        //Update tau boundary condition, added 27.08. PJS
        forAll(tau()().boundaryField(), patchI)
        {
            if
            (
                (
                    tau()().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<symmTensor>::typeName
                )
                ||
                (
                    tau()().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<symmTensor>::typeName
                )
                ||
                (
                    tau()().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<symmTensor>::typeName
                )
            )
            {
                correctedFvPatchField<symmTensor>& aTau =
                    refCast<correctedFvPatchField<symmTensor> >
                    (
                        tau()().boundaryField()[patchI]
                    );

                aTau.movePatchSubMesh();
            }
        } 

    }


    return true;
}


bool freeSurface::moveMeshPoints()
{
        scalarField sweptVolCorr =
            phi_.boundaryField()[aPatchID()]
          - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

        word ddtScheme
        (
            mesh().schemesDict().ddtScheme
            (
                "ddt(" + rho().name() + ',' + U().name()+')'
            )
        );

        if
        (
            ddtScheme
         == fv::CrankNicolsonDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else if
        (
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
        else
        {
            FatalErrorIn("freeSurface::movePoints()")
                << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
        }


        const scalarField& Sf = aMesh().S();
        const vectorField& Nf = aMesh().faceAreaNormals().internalField();

        scalarField deltaH =
            sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));


        pointField displacement = pointDisplacement(deltaH);


        //-- Set mesh motion boundary conditions

        tetPointVectorField& motionU =
            const_cast<tetPointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        fixedValueTetPolyPatchVectorField& motionUaPatch =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        tetPolyPatchInterpolation tppiAPatch
        (
            refCast<const faceTetPolyPatch>
            (
                motionUaPatch.patch()
            )
        );

        motionUaPatch ==
            tppiAPatch.pointToPointInterpolate
            (
                displacement/DB().deltaT().value()
            );

        if (twoFluids())
        {
            fixedValueTetPolyPatchVectorField& motionUbPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            tetPolyPatchInterpolation tppiBPatch
            (
                refCast<const faceTetPolyPatch>(motionUbPatch.patch())
            );

            motionUbPatch ==
                tppiBPatch.pointToPointInterpolate
                (
                    interpolatorAB().pointInterpolate
                    (
                        displacement/DB().deltaT().value()
                    )
                );
        }

        mesh().update();

        // faMesh motion is done automatically, using meshObject
        // HJ, 8/Aug/2011
//         aMesh().movePoints(mesh().points());


        // Move correctedFvPatchField fvSubMeshes                                                   //

        forAll(U().boundaryField(), patchI)
        {
            if
            (
                (
                    U().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<vector>::typeName
                )
                ||
                (
                    U().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<vector>::typeName
                )
            )
            {
                correctedFvPatchField<vector>& aU =
                    refCast<correctedFvPatchField<vector> >
                    (
                        U().boundaryField()[patchI]
                    );

                aU.movePatchSubMesh();
            }
        }

        forAll(p().boundaryField(), patchI)
        {
            if
            (
                (
                    p().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<scalar>::typeName
                )
                ||
                (
                    p().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<scalar>::typeName
                )
            )
            {
                correctedFvPatchField<scalar>& aP =
                    refCast<correctedFvPatchField<scalar> >
                    (
                        p().boundaryField()[patchI]
                    );

                aP.movePatchSubMesh();
            }
        }

        //Update tau boundary conditions, added 27.08. PJS
        forAll(tau()().boundaryField(), patchI)
        {
            if
            (
                (
                    tau()().boundaryField()[patchI].type()
                 == fixedGradientCorrectedFvPatchField<symmTensor>::typeName
                )
                ||
                (
                    tau()().boundaryField()[patchI].type()
                 == fixedValueCorrectedFvPatchField<symmTensor>::typeName
                )
                ||
                (
                    tau()().boundaryField()[patchI].type()
                 == zeroGradientCorrectedFvPatchField<symmTensor>::typeName
                )
            )
            {
                correctedFvPatchField<symmTensor>& aTau =
                    refCast<correctedFvPatchField<symmTensor> >
                    (
                        tau()().boundaryField()[patchI]
                    );

                aTau.movePatchSubMesh();
            }
        }        

    return true;
}


void freeSurface::updateBoundaryConditions()
{
    updateVelocity();
    //updateStress();
    updateSurfactantConcentration();
    updatePressure();
}

void freeSurface::updateStress()
{



    tauCorrect();

    symmTensorField nGradTau = tau()().boundaryField()[aPatchID()].snGrad();                   



    if                                                                      
    (
        tau()().boundaryField()[aPatchID()].type()
     == fixedGradientCorrectedFvPatchField<symmTensor>::typeName
    )
    {
        fixedGradientCorrectedFvPatchField<symmTensor>& aTau =
            refCast<fixedGradientCorrectedFvPatchField<symmTensor> >
            (
                tau()().boundaryField()[aPatchID()]
            );

            aTau.gradient() = nGradTau;                                 
    }
    else if
    (
        tau()().boundaryField()[aPatchID()].type()
     == fixedGradientFvPatchField<symmTensor>::typeName
    )
    {
        fixedGradientFvPatchField<symmTensor>& aTau =
            refCast<fixedGradientFvPatchField<symmTensor> >
            (
                tau()().boundaryField()[aPatchID()]
            );

        aTau.gradient() = nGradTau;
    }

}




void freeSurface::updateVelocity()
{
    if(twoFluids())
    {
        vectorField nA = mesh().boundary()[aPatchID()].nf();                    // calculate normal vector from boundary patch


        areaVectorField nAs = aMesh().faceAreaNormals();                        


        vectorField nB = mesh().boundary()[bPatchID()].nf();

        scalarField DnB = interpolatorBA().faceInterpolate                      // inverse normal distance from centroid of cell BP
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();          // gives inverse distance between cell center and boundary face

        tauCorrect(); 																// calculate tau
        tauSurfaceCorrect();

        



        


// TOTAL Patch A

        vectorField UtPA =                                                      // boundary cell centroid velocity of patch A / velocity on the cell faces v_Af
            U().boundaryField()[aPatchID()].patchInternalField();

        if                                                                      // update boundary condition
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            UtPA += aU.corrVecGrad();                                           // aU? acceleration = gradient, if there is fixedgradient
        }

// TANGENTIAL Patch A

        UtPA -= nA*(nA & UtPA);                                                 // tangential component of the velocity in the cell centroid of boundary patch A

// TOTAL Patch B

        vectorField UtPB = interpolatorBA().faceInterpolate                     // boundary velocity of patch B // velocity on cell faces
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        if                                                                      // correct boundary velocity
        (
            U().boundaryField()[bPatchID()].type()
         == fixedValueCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedValueCorrectedFvPatchField<vector>& bU =
                refCast<fixedValueCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[bPatchID()]
                );

            UtPB += interpolatorBA().faceInterpolate(bU.corrVecGrad());         // velocity of patch B interpolated?
        }

// TANGENTIAL Patch B

        UtPB -= nA*(nA & UtPB);                                                 // tangential boundary velocity = boundary velocity - normal velocity

        vectorField UtFs = muFluidA().value()*DnA*UtPA                          // why mu*Dn*Ut? the discretisation says (mu/Dn)*Ut
          + muFluidB().value()*DnB*UtPB;                                        

// NORMAL

        vectorField UnFs =                                                      // UnFs == normal free surface tension force
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();                              // magSf = magnitide of the cell surface, UnFs = n*phi/A = flux

        Us().internalField() += UnFs - nA*(nA & Us().internalField());            // Us == boundary velocity, alte normale geschwinigkeit abziehen (na*Us); neue berechnete addireren Unfs
        
        correctUsBoundaryConditions();
 
// TANGENTIAL

        UtFs -= (muFluidA().value() - muFluidB().value())*                      
            ( fac::grad( Us() ) & aMesh().faceAreaNormals() )().internalField();  


        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if(!cleanInterface())
        {
            tangentialSurfaceTensionForce =
                surfaceTensionGrad()().internalField();
        }
        else                                                                    // without surfactants
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()                             // sigma = const (0.072 Nm)
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()                 // aMesh().Le == calculate edge length; integral surface tension force over edges
                )().internalField();                                            // 

            tangentialSurfaceTensionForce =                                     // see Brackbill et al. Modeling Surface Tensiom F_sa = sigma * kappa * n
                surfaceTensionForce
              - cleanInterfaceSurfTension().value()                             // faceCurvatures in faMeshDemandDrivenData.C
               *aMesh().faceCurvatures().internalField()*nA;
        }

        UtFs += tangentialSurfaceTensionForce;                                  // total tangential force: grad. vel + surface tension; Here, UtFs = force

        UtFs -= (nAs & tauSurface());                                                   // added: viscoelastic stress term

        UtFs += (nAs * (nAs & (nAs & tauSurface())));

        UtFs /= muFluidA().value()*DnA + muFluidB().value()*DnB + VSMALL;       // VSMALL != 0; Here, UtFs = velocity. 

// NORMAL + TANGENTIAL
        Us().internalField() = UnFs + UtFs;                                     // 
        correctUsBoundaryConditions();

        // Store old-time velocity field U()
        U().oldTime();

        U().boundaryField()[bPatchID()] ==
            interpolatorAB().faceInterpolate(UtFs)
          + nB*fvc::meshPhi( rho(),U() ) ().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf();

        if
        (
            p().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& pB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    p().boundaryField()[bPatchID()]
                );

            pB.gradient() =
               - rhoFluidB().value()
                *(
                     nB&fvc::ddt(U())().boundaryField()[bPatchID()]
                 );
        }


        // Update fixedGradient boundary condition on patch A

        vectorField nGradU =
            muFluidB().value()*(UtPB - UtFs)*DnA
          + tangentialSurfaceTensionForce
          - muFluidA().value()*nA*fac::div(Us())().internalField()
          + (muFluidB().value() - muFluidA().value())
           *(fac::grad(Us())().internalField()&nA);

        nGradU -= (nAs & tauSurface());                                         // added: vicsoelastic stress term

        nGradU += (nAs * (nAs & (nAs & tauSurface())));

        nGradU /= muFluidA().value() + VSMALL; 

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name()
                    <<  " for freeSurface patch is "
                    << U().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << fixedGradientCorrectedFvPatchField<vector>::typeName
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName
                    << abort(FatalError);
        }
    }


    // -if (!= twoFluids()) do...                                   ****************************************NOT IMPORTANT FROM HERE*************************************
    else
    {
        vectorField nA = aMesh().faceAreaNormals().internalField();

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        // Correct normal component of free-surface velocity
        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if(!cleanInterface())
        {
            tangentialSurfaceTensionForce =
                surfaceTensionGrad()().internalField();
        }
        else
        {
            vectorField surfaceTensionForce =
                cleanInterfaceSurfTension().value()
               *fac::edgeIntegrate
                (
                    aMesh().Le()*aMesh().edgeLengthCorrection()
                )().internalField();

            tangentialSurfaceTensionForce =
                surfaceTensionForce
              - cleanInterfaceSurfTension().value()
               *aMesh().faceCurvatures().internalField()*nA;

            if (muFluidA().value() < SMALL)
            {
                tangentialSurfaceTensionForce = vector::zero;
            }
        }

        vectorField tnGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField UtPA =
            U().boundaryField()[aPatchID()].patchInternalField();
        UtPA -= nA*(nA & UtPA);

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();           

        vectorField UtFs = UtPA + tnGradU/DnA;

        Us().internalField() = UtFs + UnFs;
        correctUsBoundaryConditions();

        vectorField nGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - nA*fac::div(Us())().internalField()
          - (fac::grad(Us())().internalField()&nA);

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name()
                    <<  " for freeSurface patch is "
                    << U().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << fixedGradientCorrectedFvPatchField<vector>::typeName
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName
                    << abort(FatalError);
        }
    }
}


void freeSurface::updatePressure()
{
    // Correct pressure boundary condition at the free-surface

    tauCorrect(); 																			// calculate tau
    tauSurfaceCorrect();
    //symmTensorField tauPA = tau()().boundaryField()[aPatchID()].patchInternalField();   //tau boundary field of patch A

    vectorField nA = mesh().boundary()[aPatchID()].nf();                                // normal vector in free surface

    areaVectorField nAs = aMesh().faceAreaNormals();                                    // addded for fa calculations


    if(twoFluids())
    {
        scalarField pA =                                                                // pA = pB
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();               // free surface curvature 

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;

        if(cleanInterface())                                                            // clean interface?
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;                                     // surfTension.internalField = sigma; 

            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*(muFluidA().value() - muFluidB().value())                             
            *fac::div(Us())().internalField();

        pA -= (nAs & (nAs & (tauSurface())));                                                    // added.

//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
        vector R0 = vector::zero;

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        p().boundaryField()[aPatchID()] == pA;                                          
    }
    else
    {
//         vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
        vector R0 = vector::zero;

        scalarField pA =
          - rhoFluidA().value()*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K)
            << endl;

        if(cleanInterface())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            scalarField surfTensionK =
                surfaceTension().internalField()*K;

            pA -= surfTensionK - gAverage(surfTensionK);
        }

        pA -= 2.0*muFluidA().value()*fac::div(Us())().internalField();

        pA -= (nAs & (nAs & (tauSurface())));                                                    // added.

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed apsolute
    // pressure

//     vector R0 = gAverage(mesh().C().boundaryField()[aPatchID()]);
    vector R0 = vector::zero;

    for (int patchI=0; patchI < p().boundaryField().size(); patchI++)                   // update boundary condition
    {
        if
        (
            p().boundaryField()[patchI].type()
         == fixedValueFvPatchScalarField::typeName
        )
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] ==
                  - rho().boundaryField()[patchI]
                   *(g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}


void freeSurface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void freeSurface::updateSurfactantConcentration()
{
    if(!cleanInterface())
    {
        Info << "Correct surfactant concentration" << endl << flush;

        updateSurfaceFlux();

        // Crate and solve the surfactanta transport equation
        faScalarMatrix CsEqn
        (
            fam::ddt(surfactantConcentration())
          + fam::div(Phis(), surfactantConcentration())
          - fam::laplacian
            (
                surfactant().surfactDiffusion(),
                surfactantConcentration()
            )
        );


        if(surfactant().soluble())
        {
            const scalarField& C =
                mesh().boundary()[aPatchID()]
               .lookupPatchField<volScalarField, scalar>("C");

            areaScalarField Cb
            (
                IOobject
                (
                    "Cb",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                dimensioned<scalar>("Cb", dimMoles/dimVolume, 0),
                zeroGradientFaPatchScalarField::typeName
            );

            Cb.internalField() = C;
            Cb.correctBoundaryConditions();

            CsEqn +=
                fam::Sp
                (
                    surfactant().surfactAdsorptionCoeff()*Cb
                  + surfactant().surfactAdsorptionCoeff()
                   *surfactant().surfactDesorptionCoeff(),
                    surfactantConcentration()
                )
              - surfactant().surfactAdsorptionCoeff()
               *Cb*surfactant().surfactSaturatedConc();
        }

        CsEqn.solve();

        Info << "Correct surface tension" << endl;

        surfaceTension() =
            cleanInterfaceSurfTension()
          + surfactant().surfactR()
           *surfactant().surfactT()
           *surfactant().surfactSaturatedConc()
           *log(1.0 - surfactantConcentration()
           /surfactant().surfactSaturatedConc());

        if(neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn
            (
                "void freeSurface::correctSurfactantConcentration()"
            )   << "Surface tension is negative"
                    << abort(FatalError);
        }
    }
}


void freeSurface::correctUsBoundaryConditions()
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryField()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if(ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        isA<slipFvPatchVectorField>
                        (
                            U().boundaryField()[ngbPolyPatchID]
                        )
                    )
                 ||
                    (
                        isA<symmetryFvPatchVectorField>
                        (
                            U().boundaryField()[ngbPolyPatchID]
                        )
                    )
                )
                {
                    vectorField N =
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


vector freeSurface::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return gSum(pressureForces);
}


vector freeSurface::totalViscousForce() const
{
    const scalarField& S = aMesh().S();                                     // surface
    const vectorField& n = aMesh().faceAreaNormals().internalField();       // Area normals
    const areaVectorField& nAs = aMesh().faceAreaNormals();                        


    //tauSurfaceCorrect();                                                     

    vectorField nGradU =                                                    // velocity gradient on face
        U().boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces =
      - muFluidA().value()*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        );


    viscousForces += nAs & tauSurface();                                                   // added: viscoelastic stress term

    return gSum(viscousForces);
}


vector freeSurface::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();

    vectorField surfTensionForces(n.size(), vector::zero);

    if(cleanInterface())
    {
        surfTensionForces =
            S*cleanInterfaceSurfTension().value()
           *fac::edgeIntegrate
            (
                aMesh().Le()*aMesh().edgeLengthCorrection()
            )().internalField();
    }
    else
    {
        surfTensionForces *= surfaceTension().internalField()*K;
    }

    return gSum(surfTensionForces);
}


void freeSurface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();


    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll(faces, faceI)
    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), vector::zero);

    forAll (faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].normal(newPoints);
    }

    forAll(deltaH, faceI)
    {
        deltaH[faceI] = sweptVol[faceI]/
            (faceArea[faceI] & facesDisplacementDir()[faceI]);
    }

    displacement = pointDisplacement(deltaH);
}


scalar freeSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if(cleanInterface())
    {
        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/(cleanInterfaceSurfTension().value() + SMALL)
            )
        );
    }
    else
    {
        scalarField sigmaE =
            linearEdgeInterpolate(surfaceTension())().internalField()
          + SMALL;

        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}


void freeSurface::updateProperties()                                                                
{
    
    
	muFluidA_ = dimensionedScalar("muFluidA",dimensionSet(1,-1,-1,0,0,0,0), scalar(0));

    Info << "Calculating muFluidA" << endl;

    //Sum up solvent related viscosities
    PtrList<entry> modelEntries_(subDict("rheology").lookup("models"));                         // added PJS, 19.08.19

	Info << "Moldel entry list created" << endl;

    forAll(modelEntries_, modelI)
    {
        muFluidA_ += dimensionedScalar(modelEntries_[modelI].dict().lookup("etaS")); 
    }
    Info << "muFluidA is calculated as:" << muFluidA_ << endl;



	//muFluidA_ = dimensionedScalar(this->lookup("muFluidA"));

    muFluidB_ = dimensionedScalar(this->lookup("muFluidB"));

    rhoFluidA_ = dimensionedScalar(this->lookup("rhoFluidA"));

    rhoFluidB_ = dimensionedScalar(this->lookup("rhoFluidB"));

    g_ = dimensionedVector(this->lookup("g"));

    cleanInterfaceSurfTension_ =
        dimensionedScalar(this->lookup("surfaceTension"));
}


void freeSurface::writeVTK() const
{
    aMesh().patch().writeVTK
    (
        DB().timePath()/"freeSurface",
        aMesh().patch(),
        aMesh().patch().points()
    );
}


void freeSurface::writeVTKControlPoints()
{
    // Write patch and points into VTK
    fileName name(DB().timePath()/"freeSurfaceControlPoints");
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << controlPoints().size() << " float" << nl;

    forAll(controlPoints(), pointI)
    {
        mps << controlPoints()[pointI].x() << ' '
            << controlPoints()[pointI].y() << ' '
            << controlPoints()[pointI].z() << nl;
    }

    // Write vertices
    mps << "VERTICES " << controlPoints().size() << ' '
        << controlPoints().size()*2 << nl;

    forAll(controlPoints(), pointI)
    {
        mps << 1 << ' ' << pointI << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

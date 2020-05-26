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

Application
    interfaceTrackinFoam

Description
    Incompressible laminar CFD code for simulation of a single bubble rising
    in a stil liquid. Interface between fluid phases is tracked using moving
    mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "dynamicFvMesh.H"
#include "viscoelasticLaw.H"                                                                        // added
#include "autoPtr.H"          
#include "freeSurface.H"
#include "inletOutletFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "inletOutletFaPatchFields.H"
#include "fixedGradientFaPatchFields.H"
#include "boundBox.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    pimpleControl pimple(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createBubble.H"
#   include "createSurfactantConcentrationField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

        interface.moveMeshPointsForOldFreeSurfDisplacement();               // freeSurface Method: bool method checks for mesh, fluid...; object created in CreateFields from Class FreeSurface defined in FreeSurface.C

        interface.updateDisplacementDirections();                           // freeSurface Method: void method         
                                                                            // Correcte point displacement direction
                                                                            // at the "centerline" symmetryPlane which represents the axis
                                                                            // of an axisymmetric case

        interface.predictPoints();                                          // freeSurface Method: bool method, freeSurfCorr-loop to smooth interface 

        Info<< "\nMax surface Courant Number = "
            << interface.maxCourantNumber() << endl << endl;

        while (pimple.loop())
        {
            // Update interface bc
            interface.updateBoundaryConditions();                           //freeSurface Method: updates velocity, surfactantConcentration, pressure

            Info<< "boundary conditions updated" << endl << endl;


            // Make the fluxes relative
            phi -= fvc::meshPhi(rho, U);                                    // phi: interpolated velocity normal to cell surface; calculates new phi

#           include "CourantNo.H"                                           // calculating courant no: is a necessary condition for convergence while solving certain partial differential equations (usually hyperbolic PDEs) numerically.                                      

             interface.tauCorrect();                                            // added;

             Info<< "tau corrected" << endl << endl;


#           include "UEqn.H" 

            UEqn.relax();                                               // relatxation for momentum eqn

            solve(UEqn == -fvc::grad(p));

            Info<< "UEqn solved" << endl << endl;

            // --- PISO loop
            while (pimple.correct())
            {
                volScalarField AU = UEqn.A();

                U = UEqn.H()/AU;

                phi = (fvc::interpolate(U) & mesh.Sf());

#               include "scalePhi.H"

#               include "scaleSpacePhi.H"


                // store previous pressure for unter-relaxation, PJS added 31.07.19
                p.storePrevIter();

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
#                   include "pEqn.H"

#                   include "setReference.H"

                    pEqn.relax();                           // PJS, added 30.07.19

                    pEqn.solve();

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"


                // Explicitly relax pressure for momentum corrector
                p.relax();
                
                // Momentum corrector
                U -= fvc::grad(p)/AU;
                U.correctBoundaryConditions();
            }

#           include "solveBulkSurfactant.H"

            interface.correctPoints();

#           include "freeSurfaceContinuityErrs.H"
        }

#       include "updateMovingReferenceFrame.H"                      // calculates bubble centre, velocity, acceleration and position from mesh.V() and mesh.C() 

#       include "volContinuity.H"

        Info << "Total surface tension force: "
            << interface.totalSurfaceTensionForce() << endl;

        vector totalForce =                                         // summing up forces
            interface.totalViscousForce()
          + interface.totalPressureForce();

        Info << "Total force: " << totalForce << endl;


        tauTotal = interface.calculateTauTotalField(U, fluidIndicator);        // added, PJS 07.08.19: caluclate total tau field for write command

        runTime.write();

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //

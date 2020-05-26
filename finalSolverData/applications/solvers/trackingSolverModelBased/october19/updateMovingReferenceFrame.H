{
    vector Cb0 = Cb;                                                        // Cb == bubble centre

    if (interface.twoFluids())
    {
        Cb = gSum((1.0 - fluidIndicator.internalField())*                   // gSum == sums over all processors in a parallel run; fluidIndicator.internalField is 1 or 0
            mesh.C().internalField()*mesh.V())/                             // depending on liquid or bubble; 
            gSum((1.0 - fluidIndicator.internalField())*mesh.V());          // mesh.V(), mesh.C() ???
    }
    else
    {
        scalar Va = gSum(mesh.V().field());
        vector Ra = gSum(mesh.C().internalField()*mesh.V().field())/Va;

        scalar Vb =
          - gSum
            (
                mesh.Cf().boundaryField()[interface.aPatchID()]
              & mesh.Sf().boundaryField()[interface.aPatchID()]
            )/3;

        scalar V = Va + Vb;

        Cb = (V*Cbf - Va*Ra)/Vb;

        if (mesh.nGeometricD() != 3)
        {
            FatalErrorIn(args.executable())
                << "One-fluid bubble centroid calc "
                    << "is not implemented for 2d mesh"
                    << abort(FatalError);
        }
    }

    Info << "\nCurrent bubble centre: " << Cb << endl;

    dimensionedVector dUF                                               // bubble velocity difference
    (
        "dUF",
        dimVelocity,
        lambdaFf*(Cb-Cbf)/runTime.deltaT0().value()                     // lambdaFf  && lambdaF0 == 0, see constant/mrfProperties
      + lambdaF0*(Cb-Cb0)/runTime.deltaT0().value()
    );

    XF += (UF + 0.5*dUF)*runTime.deltaT0();                             // change in postion due to velocity*time, why 0.5 velociy
    Info << "Current bubble position: " << XF << endl;

//     dimensionedVector UF00 ("UF00", dimVelocity, UF0.value());
//     UF0 = UF;
    UF += dUF;                                                          // new velocity
    Info << "Current bubble velocity: " << UF << endl;

    aF = dUF/runTime.deltaT0();                                         // bubble acceleration as a result of velocity/dt
//     aF = (3.0*UF - 4.0*UF0 + UF00)/(2.0*runTime.deltaT0());
    Info << "Current bubble acceleration: " << aF << endl;

    if
    (
        U.boundaryField()[spacePatchID].type()
     == inletOutletFvPatchVectorField::typeName
    )
    {
        inletOutletFvPatchVectorField& spaceU =
            refCast<inletOutletFvPatchVectorField>
            (
                U.boundaryField()[spacePatchID]
            );

        spaceU.refValue() = -UF.value();
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Velocity boundary condition at space patch is not "
                << inletOutletFvPatchVectorField::typeName
                << abort(FatalError);
    }

    // Update movingReferenceFrame dictionary
    //- bubble centre
    if (movingReferenceFrame.found("Cb"))
    {
        movingReferenceFrame.remove("Cb");
    }
    movingReferenceFrame.add("Cb", Cb);

    //- bubble centre at the start of the caluclation
    if (movingReferenceFrame.found("Cbf"))
    {
        movingReferenceFrame.remove("Cbf");
    }
    movingReferenceFrame.add("Cbf", Cbf);

    //- bubble position vector
    if (movingReferenceFrame.found("XF"))
    {
        movingReferenceFrame.remove("XF");
    }
    movingReferenceFrame.add("XF", XF);

    //- bubble velocity
    if (movingReferenceFrame.found("UF")) 
    {
        movingReferenceFrame.remove("UF");
    }
    movingReferenceFrame.add("UF", UF);

    //- bubble acceleration
    if (movingReferenceFrame.found("aF"))
    {
        movingReferenceFrame.remove("aF");
    }
    aF.name() = "aF";
    movingReferenceFrame.add("aF", aF);
}

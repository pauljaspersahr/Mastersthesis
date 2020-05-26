fvVectorMatrix UEqn                                       
(
fvm::ddt(rho, U)
+ rho*aF
+ fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")



- fvm::laplacian(interface.muFluidB(), U)															// tau fluid B = inside
- (fluidIndicator * (interface.divTau(U) - fvm::laplacian(interface.muFluidB(), U))) 				// tau fluid A = outside



);





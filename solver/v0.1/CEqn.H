{
    volScalarField& diffEff = diffCo;
    volScalarField& alpha = alpha_solid;
    
    forAll(mesh.C(),i)
    {
	if(alpha[i]>0.5)
	{
		diffEff[i] = 0.0;
	}
    }
 
    Info<< "min/max(diffEff) = " << min(diffEff).value() << ", " << max(diffEff).value() << endl;


    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi, C, "div(phi,T)")
      - fvm::laplacian(diffEff, C,  "laplacian(kappa,T)")
      ==
       // fluid.heatTransfer(C)
      fvOptions(C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);

    CEqn.solve();

    fvOptions.correct(C);

    Info<< "min/max(C) = "
        << min(C).value() << ", " << max(C).value() << endl;
}

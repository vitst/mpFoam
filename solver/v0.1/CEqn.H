{
    volScalarField& diffEff = diffCo;
    phaseModel& solidPhase = fluid.phase(1);
    const volScalarField& alpha_solid = solidPhase;
    //mesh.lookupObject<volScalarField>("alpha.solid").oldTime();
    
    forAll(mesh.C(),i)
    {
	if(alpha_solid[i]>0.5)
	{
		diffEff[i] = 0.0;
	}
    }
 
    Info<< "min/max(diffEff) = " << min(diffEff).value() << ", " << max(diffEff).value() << endl;

//    Su = fluid.preRect(C,molW) & C;

    forAll(mesh.C(),i)
    {
	bool flag = true;
	forAll(mesh.cellCells()[i], j)
	{
		if(alpha_solid[mesh.cellCells()[i][j]]>0.9)
		{
			flag = false;
			break;
		}
	}
	if(!flag)
	{
		diffEff[i] = 0.0;
	}
    }

//    Info<< "min/max(Su) = " << min(Su).value() << ", " << max(Su).value() << endl;


    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi, C, "div(phi,T)")
      - fvm::laplacian(diffEff, C,  "laplacian(kappa,T)")
      ==
        fluid.preRect(C,molW,SuOut) 
      + fvOptions(C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);

    CEqn.solve();

    fvOptions.correct(C);

    Info<< "min/max(C) = "
        << min(C).value() << ", " << max(C).value() << endl;
}

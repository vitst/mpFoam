{
    Info<< "Compute concentration field ..." << endl;
    volScalarField& diffEff = diffCo;
    const volScalarField& liquidPhase = alpha1;

    //- update diffusion coefficient basing on concentration and ion cloud
    
    if(ionEffect)
    {
        Info<< "Correct diffusion coefficients for ionic strength ..." << endl;
        forAll(mesh.C(),i)
        {
            if(C[i]<=0.0)
            {
                diffEff[i] = 0;
            }
            else
            {
                scalar CStar = Foam::sqrt(C[i]/conSat);
                diffEff[i] = diffInLiquid.value()*(1.0-(ionEffectA*CStar/2.0/pow((1+ionEffectB*CStar),2)))*Foam::pow(liquidPhase[i],2);
            }
            if(diffEff[i]<1e-10)
            {
                diffEff[i]=0.0;
            }
        }
    }
    else
    {
        forAll(mesh.C(),i)
        {
            diffEff[i] = diffInLiquid.value()*pow(liquidPhase[i],2);
            if(diffEff[i]<1e-10)
            {
                diffEff[i]=0.0;
            }
        }
    }
    
 
    Info<< "min/max(diffEff) = " << min(diffEff).value() << ", " << max(diffEff).value() << endl;

//    Su = fluid.preRect(C,molW) & C;

//    forAll(mesh.C(),i)
//    {
//	bool flag = true;
//	forAll(mesh.cellCells()[i], j)
//	{
//		if(alpha_solid[mesh.cellCells()[i][j]]>0.9)
//		{
//			flag = false;
//			break;
//		}
//	}
//	if(!flag)
//	{
//		diffEff[i] = 0.0;
//	}
//    }

//    Info<< "min/max(Su) = " << min(Su).value() << ", " << max(Su).value() << endl;


    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi, C, "div(phi,T)")
      - fvm::laplacian(diffEff, C,  "laplacian(kappa,T)")
      ==
        phaseChange.preRect(C, SuOut)
      //- phaseChange.Kexp(C)/molW
      + fvOptions(C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);

    CEqn.solve();

    fvOptions.correct(C);

    forAll(mesh.C(),i)
    {
        if(C[i]<0.0)
        {
            C[i]=0.0;
        }
    }

    //fluid.correct();

    Info<< "min/max(C) = "
        << min(C).value() << ", " << max(C).value() << endl;

}

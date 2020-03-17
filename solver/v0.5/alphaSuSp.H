volScalarField Su = - phaseChange.Kexp(C).ref()/rho1;
Su += fvc::div(phi)*min(max(alpha1, scalar(0)), scalar(1));
Info << "Su max/min: " << max(Su).value() << ", " << min(Su).value() << endl;
//zeroField Su;

zeroField Sp;

zeroField divU;

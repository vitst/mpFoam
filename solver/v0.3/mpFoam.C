/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mpFoam

Group
    grpMultiphaseSolvers

Description
    Solver for N incompressible, non-isothermal immiscible fluids with
    phase-change.  Uses a VOF (volume of fluid) phase-fraction based interface
    capturing approach.

    The momentum, energy and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "dynamicFvMesh.H"
#include "subCycle.H"
#include "multiphaseSystem.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"
//#include "radiationModel.H"
#include "HashPtrTable.H"
#include "fvcDDt.H"
#include "zeroField.H"

//#include "CorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for N incompressible, non-isothermal immiscible fluids with"
        " phase-change,"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    //- Added for mesh refinement
    //#include "createDynamicFvMesh.H"
    //#include "initContinuityErrs.H"
    //#include "createDyMControls.H"
    
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1.0)
    );

    #include "createFieldRefs.H"
    //#include "initCorrectPhi.H"
    //#include "createUf.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //#include "readDyMControls.H"
        //volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        fluid.solve();
        rho = fluid.rho();

        Info<< "Start pimple loop ... " << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            
            solve(fvm::ddt(rho) + fvc::div(rhoPhi));
            #include "UEqn.H"
            //#include "YEqns.H"
            //#include "TEqn.H"
	        #include "CEqn.H"
            fluid.correct();
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            //Info<< "pimple correct finished" << endl;
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        
        Info<< "pimp loop finished!" << endl;

        rho = fluid.rho();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

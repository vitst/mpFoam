#include "createRDeltaT.H"

Info<< "Reading personlized input parameters" << endl;
IOdictionary inputParameters
(
    IOobject
    (
        "inputParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar molW("Mv", dimMass/dimMoles, inputParameters);
dimensionedScalar Csat("Cactivate", dimMoles/dimVolume, inputParameters);
bool nucleationFlag(inputParameters.lookupOrDefault<bool>("nucleation", false));

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField SuOut
(
    IOobject
    (
        "SuOut",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Diffusion coefficient of acid in fluid phase
dimensionedScalar diffInLiquid
(
    "diffInLiquid",
    dimLength*dimLength/dimTime,
    inputParameters
);


Info<< "Reading field diffusion coefficient for C field\n" << endl;
volScalarField diffCo
(
    IOobject
    (
        "diffCo",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field for nucleation\n" << endl;
volScalarField nuSite
(
    IOobject
    (
        "nuSite",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

volScalarField nuRate
(
    IOobject
    (
        "nuRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimArea/dimTime, Zero)
);

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

dimensionedScalar liquidMass_(0.0);
dimensionedScalar liquidMassOld_(0.0);
dimensionedScalar solidMass_(0.0);
dimensionedScalar dmdtTot_(0.0);
dimensionedScalar deltaT_(0.0);

vectorList nuSiteList;
dimensionedScalar nuCount_(0.0);
dimensionedScalar nuTotal_(0.0);

scalar timeMarker = 0.0;

volScalarField Cmask
(
    IOobject
    (
        "Cmask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*0.0
);

phaseChangeReaction phaseChange = phaseChangeReaction(inputParameters, mesh, C, alpha1, Cmask);

//- diffusion coefficient flag and inputs
bool ionEffect = inputParameters.lookupOrDefault<bool>("ionEffect", false);
scalar ionEffectA = inputParameters.lookupOrDefault<scalar>("ionEffectA", 0.0);
scalar ionEffectB = inputParameters.lookupOrDefault<scalar>("ionEffectB", 0.0);
scalar conSat = inputParameters.lookupOrDefault<scalar>("Cactivate", 1.0);

//- initial nucleation vof
scalar nuIntVof = inputParameters.lookupOrDefault<scalar>("nucleationVof", 1e-6);
scalar nuIntVofCut = inputParameters.lookupOrDefault<scalar>("nucleationVofCut", 0.5);

#include "createMRF.H"
#include "createFvOptions.H"

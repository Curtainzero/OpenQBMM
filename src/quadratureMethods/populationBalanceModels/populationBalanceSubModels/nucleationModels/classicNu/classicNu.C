/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "classicNu.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

defineTypeNameAndDebug(classicNu, 1);

addToRunTimeSelectionTable
(
    nucleationModel,
    classicNu,
    dictionary
);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::populationBalanceSubModels::nucleationModels::classicNu::classicNu
(
    const dictionary& dict,
    const fvMesh& mesh
):
    nucleationModel(dict, mesh),
    A_("parameterA",inv(dimVolume*dimTime),dict),
    B_("parameterB",dimless,dict),
    moleDensity_("moleDensity",dimMass*inv(dimMoles),dict),
    particleRho_("rho",inv(dimVolume)*dimMass,dict)
{
    Info << A_ <<endl;
    Info << B_ <<endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::populationBalanceSubModels::nucleationModels::classicNu::~classicNu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar 
Foam::populationBalanceSubModels::nucleationModels::classicNu::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const {
    const volScalarField& totalN_ = mesh_.lookupObject<volScalarField>("moment.0.populationBalance");
    const volScalarField& SuperS_ = mesh_.lookupObject<volScalarField>("SuperS");
    const volScalarField& dcdt_ = mesh_.lookupObject<volScalarField>("dcdt");
    //Info << SuperS_[celli] <<endl;
    scalar A = A_.value();
    scalar B = B_.value();
    Foam::scalar tmp = exp(-1*(B/sqr(log(SuperS_[celli]+1e-10))));
    Foam::scalar result = A*SuperS_[celli]*tmp;
    //Info << result <<endl;
    scalar minSize = 5e-9;

    scalar kv(Foam::constant::mathematical::pi/6);

    if (((dcdt_[celli]*moleDensity_/particleRho_).value()-mesh_.V()[celli]*result*kv*pow(minSize,3)) <0) {
        result = (dcdt_[celli]*moleDensity_/particleRho_).value() / (mesh_.V()[celli]*kv*pow(minSize,3));
    }

    //if (totalN_[celli]> 100) {
        //result = 0.0;
    //}


    return (
        result*pow(minSize, momentOrder)
    );
}


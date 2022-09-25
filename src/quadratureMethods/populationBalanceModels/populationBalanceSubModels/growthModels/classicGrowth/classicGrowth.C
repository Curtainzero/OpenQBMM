/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi and 2017 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
-------------------------------------------------------------------------------
2017-03-28 Alberto Passalacqua: Adapted to single scalar calculation.
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

#include "classicGrowth.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "nucleationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(classicGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        classicGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::classicGrowth
::classicGrowth
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    growthModel(dict, mesh),
    moleDensity_("moleDensity",dimMass*inv(dimMoles),dict),
    particleRho_("rho",inv(dimVolume)*dimMass,dict),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", scalar(0))),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", GREAT)),
    dcdt_(mesh.lookupObject<volScalarField>("dcdt")),
    totalN_(mesh.lookupObject<volScalarField>("moment.0.populationBalance"))
{
    nucleationModel_ =
            Foam::populationBalanceSubModels::nucleationModel::New
            (
                dict.subDict("nucleationModel"),
                mesh
            );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::classicGrowth
::~classicGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// not used
Foam::scalar
Foam::populationBalanceSubModels::growthModels::classicGrowth::Kg
(
    const scalar& abscissa,
    const bool lengthBased,
    const label environment
) const
{
    
    return 0.0;
}

//重新定义 growth kernel kg函数  

Foam::scalar
Foam::populationBalanceSubModels::growthModels::classicGrowth::Kg_new
(
    const scalar& abscissa,
    const label celli,
    const bool lengthBased,
    const label environment
) {
    Foam::scalar kv(3.141592653/6);
    Foam::scalar tmp = pos0((dcdt_[celli]*moleDensity_/particleRho_).value()-(mesh_.V()[celli]*(nucleationModel_->nucleationSource(1, celli))*kv));
    Foam::scalar result = (tmp/(3*kv*totalN_[celli]));
    return result; 
}


//覆盖掉父类中phaseSpaceConvection
Foam::scalar
Foam::populationBalanceSubModels::growthModels::classicGrowth::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    scalar gSource(0);

    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return gSource;
    }
    label sizeOrder = momentOrder[sizeIndex];


    bool lengthBased = nodes[0].lengthBased();
    bool volumeFraction = nodes[0].useVolumeFraction();
    if (volumeFraction)
    {
        if (lengthBased)
        {
            sizeOrder += 3;
        }
        else
        {
            sizeOrder += 1;
        }
    }

    if (sizeOrder < 1)
    {
        return gSource;
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volScalarNode& node = nodes[pNodeI];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));
            scalar d = node.d(celli, bAbscissa);
            scalar n =
                node.n(celli, node.primaryWeight()[celli], bAbscissa);

            scalar gSourcei =
                n
               *Kg_new(d,celli, lengthBased)
               *sizeOrder
               *pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    gSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }
            gSource += gSourcei;
        }

        return gSource;
    }

    forAll(nodes, pNodeI)
    {
        const volScalarNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights()[sizeIndex], sNodei)
        {
            scalar bAbscissa =
                max
                (
                    node.secondaryAbscissae()[sizeIndex][sNodei][celli],
                    scalar(0)
                );
            scalar d = node.d(celli, bAbscissa);
            scalar n =
                node.n(celli, node.primaryWeight()[celli], bAbscissa)
               *node.secondaryWeights()[sizeIndex][sNodei][celli];

            scalar gSourcei =
                n
               *Kg_new(d,celli, lengthBased)
               *sizeOrder
               *pow(bAbscissa, sizeOrder - 1);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    gSourcei *=
                        node.secondaryWeights()[cmpt][sNodei][celli]
                       *pow
                        (
                            node.secondaryAbscissae()[cmpt][sNodei][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }
            gSource += gSourcei;
        }
    }

    return gSource;
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::classicGrowth::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    scalar gSource(0);

    const PtrList<volVelocityNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        return gSource;
    }
    label sizeOrder = momentOrder[sizeIndex];

    bool lengthBased = nodes[0].lengthBased();
    bool volumeFraction = nodes[0].useVolumeFraction();
    if (volumeFraction)
    {
        if (lengthBased)
        {
            sizeOrder += 3;
        }
        else
        {
            sizeOrder += 1;
        }
    }

    if (sizeOrder < 1)
    {
        return gSource;
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    forAll(nodes, pNodeI)
    {
        const volVelocityNode& node = nodes[pNodeI];

        scalar bAbscissa =
            max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));
        scalar d = node.d(celli, bAbscissa);
        scalar n =
            node.n(celli, node.primaryWeight()[celli], bAbscissa);

        scalar gSourcei =
            n*Kg_new(d,celli, lengthBased)*sizeOrder*pow(bAbscissa, sizeOrder - 1);

        forAll(scalarIndexes, nodei)
        {
            if (scalarIndexes[nodei] != sizeIndex)
            {
                gSourcei *=
                    pow
                    (
                        node.primaryAbscissae()[nodei][celli],
                        momentOrder[scalarIndexes[nodei]]
                    );
            }
        }
        forAll(velocityIndexes, cmpt)
        {
            gSourcei *=
                pow
                (
                    node.velocityAbscissae()[celli][cmpt],
                    momentOrder[velocityIndexes[cmpt]]
                );
        }
        gSource += gSourcei;
    }

    return gSource;
}



// ************************************************************************* //



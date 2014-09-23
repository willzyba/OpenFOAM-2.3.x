/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "solidParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<solidParticle>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        if (debug)
        {
            Info<< "Time = " << mesh_.time().timeName()
                << " trackTime = " << trackTime
                << " tEnd = " << tEnd
                << " steptFraction() = " << stepFraction() << endl;
        }

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // remember which cell the parcel is in
        // since this will change if a face is hit
        label cellI = cell();

        // trackToFace, tracks a position and returns 1.0 if the trajectory is completed
        // without hitting a face. Otherwise it stops and returns the proportion of the trajectory
        // that was completed.
        dt *= trackToFace(position() + dt*U_, td);
        // TODO: Should it be possible for dt==0? Can we react smarter when this happens.

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;

        cellPointWeight cpw(mesh_, position(), cellI, face());
        vector Uc = td.UInterp().interpolate(cpw);

        // This switch allows a zero-size particle to be defined that merely acts as a tracer
        // particle. This is usefull in combination with other flow models (ie. isoFoam) where
        // you want to track unsteady fluid motion. 
        if(d_>0)
        {        
            scalar rhoc = td.rhoInterp().interpolate(cpw);
            scalar nuc = td.nuInterp().interpolate(cpw);

            scalar rhop = td.cloud().rhop();
            scalar magUr = mag(Uc - U_);

            scalar ReFunc = 1.0;
            scalar Re = magUr*d_/nuc;

            if (Re > 0.01)
            {
                ReFunc += 0.15*pow(Re, 0.687);
            }

            // Compute the drag coefficient. There are various formulations, 
            // but this is 'similar' to a coefficient from "Clift et al, 1978", where 
            //   Cd = 24 * (1+0.15 Re^0.687) / Re           
            scalar Dc = (24.0*nuc/d_)*ReFunc*(3.0/4.0)*(rhoc/(d_*rhop));

            // Compute new velocity of particle. Initial velocity + acceleration from drag
            // This is a rearrangement to get U1, from the original equation
            //   U1 = U0 + dt*Dc*(Uc-U1) + dt*(1-rhoc/rhop)*g
            U_ = (U_ + dt*(Dc*Uc + (1.0 - rhoc/rhop)*td.g()))/(1.0 + dt*Dc);
        }
        else
        {
            // Here we are using the particles to track the fluid flow,
            // so the particle needs to move at whatever speed the underlying fluid moves at
            U_ = Uc;
        }

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                // The particle has reached a processor boundary, so we need move it to parallel processor
                td.switchProcessor = true;
            }
            else
            {
                // The particle has reached the edge of the domain, so it must vanish
                td.keepParticle = false;
            }
        } 
    }

    return td.keepParticle;
}


bool Foam::solidParticle::hitPatch
(
    const polyPatch&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::solidParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::solidParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    vector nw = tetIs.faceTri(mesh_).normal();
    nw /= mag(nw);

    scalar Un = U_ & nw;
    vector Ut = U_ - Un*nw;

    if (Un > 0)
    {
        U_ -= (1.0 + td.cloud().e())*Un*nw;
    }

    U_ -= td.cloud().mu()*Ut;
}


void Foam::solidParticle::hitPatch
(
    const polyPatch&,
    trackingData& td
)
{
    td.keepParticle = false;
}


void Foam::solidParticle::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::solidParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


Foam::scalar Foam::solidParticle::wallImpactDistance(const vector&) const
{
    if(d_<=0) return 0.0005;  // some small number. TODO: Use a beter globally defined SMALL value. 
    return 0.5*d_;
}


// ************************************************************************* //

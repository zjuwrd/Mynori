#pragma once

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
NORI_NAMESPACE_BEGIN

struct Photon
{
    public:
        Photon() {}
        explicit Photon(const Vector3f& position, const Color3f& flux, const Vector3f& dir)
        :position(position),flux(flux) {}
        
        ~Photon()=default;
        
        //location of photon
        Vector3f position=Vector3f(0.f,0.f,0.f);
        //flux carried by photon
        Color3f flux = 0.f;
        //theta phi represent the local in-direction of photon
        Vector3f dir=Vector3f(0.f,0.f,0.f);
        
};

NORI_NAMESPACE_END
#pragma once
#include<nori/object.h>
NORI_NAMESPACE_BEGIN
struct Intersection;

class Medium : public NoriObject
{
public:
    /* Compute the beam transmittance */
    virtual Color3f transmitance(const Ray3f ray, Sampler *sampler) const = 0;
    /* Sample an intersection in the medium */
    virtual Color3f sample(const Ray3f ray, Sampler *sampler, Intersection &mediumIts, Vector3f& newdir) const = 0;
    /* Wrapper of the phase function */
    virtual float SamplePhase(const Vector3f &wi, Vector3f &wo, const Point2f &sample) const = 0;
    virtual float EvalPhase(const Vector3f &wi, const Vector3f &wo) const = 0;
    virtual Color3f albedo() const = 0;
    virtual Color3f get_sigma_s() const = 0;
    
    EClassType getClassType() const { return EMedium; }


};

NORI_NAMESPACE_END

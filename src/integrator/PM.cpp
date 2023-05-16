#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/sampler.h>
#include<nori/warp.h>

NORI_NAMESPACE_BEGIN

class PM:public Integrator
{
    public:

    PM(const PropertyList &props)
    {
    
    
    }


    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        return 0.f;
    }

    std::string toString() const {
        return "PMIntegrator[]";
    }

    private:
        float radius;

};




NORI_REGISTER_CLASS(PM, "PhotonMapping");


NORI_NAMESPACE_END

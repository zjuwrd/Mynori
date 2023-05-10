#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/sampler.h>
#include<nori/warp.h>

NORI_NAMESPACE_BEGIN

class AoIntegrator:public Integrator
{
    public:

    AoIntegrator(const PropertyList &props)
    {
    
    
    }


    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        constexpr float eps = 1e-4f;

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Point2f Uniform2Dpt = sampler->next2D();
        const Vector3f wo_local = Warp::squareToCosineHemisphere(Uniform2Dpt);
        const Vector3f wo_world = its.shFrame.toWorld(wo_local);
        
        Ray3f shadowray(its.p + wo_world*eps, wo_world);

        if(scene->rayIntersect(shadowray))
            return Color3f(0.0f);
        else
            return Color3f(1.0f);

    }

    std::string toString() const {
        return "AoIntegrator[]";
    }
};




NORI_REGISTER_CLASS(AoIntegrator, "ao");


NORI_NAMESPACE_END

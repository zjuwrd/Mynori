#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
                                                                /*get the position of point light and the light intensity*/
    SimpleIntegrator(const PropertyList &props)
    :m_lightpos(props.getPoint("position")),
     m_lightIntensity(props.getColor("energy")) 
    {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        constexpr float epsilon=1e-4f;
        
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        const Vector3f dir = (m_lightpos - its.p).normalized();
        const Vector3f origin = its.p+dir*epsilon;
        Ray3f shadowRay(origin, dir);
        if(scene->rayIntersect(shadowRay))
            return Color3f(0.0f);

        auto cosTheta = its.shFrame.cosTheta(its.shFrame.toLocal(dir));

        return m_lightIntensity* std::max(0.f,cosTheta)*INV_FOURPI*INV_PI/std::pow((m_lightpos-its.p).norm(),2);
        
    }




    std::string toString() const {
        return "SimpleIntegrator[]";
    }
private:
    Point3f m_lightpos;
    Color3f m_lightIntensity;

};

NORI_REGISTER_CLASS(SimpleIntegrator,"simple");
NORI_NAMESPACE_END
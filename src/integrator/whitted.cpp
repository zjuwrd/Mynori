#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
NORI_NAMESPACE_BEGIN
class WhittedIntegrator : public Integrator {
    public:
        WhittedIntegrator(const PropertyList& props){}

#if 1
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        Intersection its;
        Color3f L(0.0f);
        if (!scene->rayIntersect(ray, its)) {
        return L;
        }
        
        //Le
        Color3f Le(0.0f);
        if (its.mesh->isEmitter()) {
        EmitterQueryRecord lRecE(ray.o, its.p, its.shFrame.n);
        Le = its.mesh->getEmitter()->eval(lRecE);
        }

        //sample direct light
        if (its.mesh->getBSDF()->isDiffuse()) {
            const Mesh* lightMesh = scene->SampleLight(sampler->next1D());
            EmitterQueryRecord lRec(its.p);
            Color3f Li = lightMesh->getEmitter()->sample( lRec, sampler);
            if (scene->rayIntersect(lRec.shadowRay())) {
                Li = 0;
            }

            float cosTheta =  std::max(0.f, Frame::cosTheta(its.shFrame.toLocal(lRec.wi)));
            BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lRec.wi), ESolidAngle);
            Color3f f = its.mesh->getBSDF()->eval(bRec);
            
            return Le + Li * f * cosTheta * float(scene->getEmissiveMeshesCount());
        } 
        else {
            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            
            Color3f refColor = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                constexpr float RR_prob=0.95f;
            if (sampler->next1D() < 0.95 && refColor.x() > 0.f) {
                return this->Li(scene, sampler, Ray3f(its.p, its.toWorld(bRec.wo))) / RR_prob * refColor;
            } 
            else 
                return 0.f;
        }
    }

#endif



        std::string toString() const {
        return "WhittedIntegrator[]";
    }


};


NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
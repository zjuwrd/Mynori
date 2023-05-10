#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
NORI_NAMESPACE_BEGIN
class WhittedIntegrator : public Integrator {
    public:
        WhittedIntegrator(const PropertyList& props){}
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            /* Find the surface that is visible in the requested direction */
            Intersection its;
            if (!scene->rayIntersect(ray, its))
                return Color3f(0.0f);

            
            Color3f Le(0.0f),Li(0.0f);
            if(its.mesh->isEmitter())
            {
                LightQueryRecord record(its.p,ray.o,its.shFrame.n);
                Le=its.mesh->getEmitter()->eval(its.mesh,record);
            }
            
            auto& meshes = scene->getMeshes();
            size_t Emmitter_count = 0;
            for(const auto& mesh:meshes)
            {
                if(mesh && mesh->isEmitter())
                {
                    Emmitter_count++;
                }
            }

            
            if(Emmitter_count<=0)
            {
                return Le;
            }
            else
            {
                for(auto mesh:meshes)
                {
                    if(mesh && mesh->isEmitter())
                    {
                        
                        LightQueryRecord record(its.p,ray.o,its.shFrame.n);
                        float pdf=0.f;
                        Color3f radiance=mesh->getEmitter()->sample(mesh,record,sampler,pdf);
                        if(pdf>0.f)
                        {
                            Vector3f wi = record.Light_Sample_point - its.p;
                            Ray3f shadow_ray(its.p, wi.normalized(), Epsilon, wi.norm() - Epsilon);
                            wi.normalize();
                            Intersection its_shadow;
                            if (!scene->rayIntersect(shadow_ray, its_shadow))
                            {
                                BSDFQueryRecord bRec(its.shFrame.toLocal(wi),its.shFrame.toLocal(-ray.d),  EMeasure::ESolidAngle);
                                Li += radiance * its.mesh->getBSDF()->eval(bRec) * std::abs(its.shFrame.n.dot(wi)) * std::abs(record.AreaLight_normal.dot(-wi))
                                /(record.Light_Sample_point - its.p).squaredNorm()/pdf;
                            }
                        }
                        
                    }
                }
                
            }

            return Le+Li;

        }

        std::string toString() const {
        return "WhittedIntegrator[]";
    }


};


NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
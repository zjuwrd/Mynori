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
                EmitterQueryRecord record(ray.o,its.p,its.shFrame.n);
                Le=its.mesh->getEmitter()->eval(record);
            }
            
            if(its.mesh->getBSDF()->isDiffuse()) // diffuse surface
            { 
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
                            
                            EmitterQueryRecord record(its.p);
                            
                            float pdf=0.f;
                            Color3f radiance=mesh->getEmitter()->sample(record,sampler);
                            pdf = record.pdf;
                            if(pdf>0.f)
                            {
                                Vector3f wi = record.wi;
                                Ray3f shadow_ray(its.p, wi, Epsilon, record.dist - Epsilon);
                                Intersection its_shadow;
                                if (!scene->rayIntersect(shadow_ray, its_shadow))
                                {
                                    BSDFQueryRecord bRec(its.shFrame.toLocal(wi),its.shFrame.toLocal(-ray.d),  EMeasure::ESolidAngle);
                                    Li += radiance * its.mesh->getBSDF()->eval(bRec) * std::abs(its.shFrame.n.dot(wi));
                                }
                            }
                            
                        }
                    }
                    
                }

                return Le+Li;
            }
            else{ //specular surface
                const auto bsdf = its.mesh->getBSDF();
                auto sample_pt = sampler->next2D();
                BSDFQueryRecord brec(its.toLocal(-ray.d));
                auto atteunance = bsdf->sample(brec,sample_pt);

                //use russian roullete
                constexpr float russian_roulete=0.95f;
                if(sampler->next1D()<russian_roulete && atteunance[0]>0.f)
                {
                    return this->Li(scene,sampler,Ray3f(its.p,its.toWorld(brec.wo))) /0.95f*atteunance ;

                }                
                else{
                    return Color3f(0.f);
                }

            }
        }

        std::string toString() const {
        return "WhittedIntegrator[]";
    }


};


NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
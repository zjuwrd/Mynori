#include<nori/integrator.h>
#include<nori/scene.h>
#include<nori/emitter.h>
#include<nori/sampler.h>
#include<nori/warp.h>
#include<nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathEms : public Integrator {
    public:
    PathEms(const PropertyList &props) {}

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_)const{
    #if 1
        constexpr uint32_t MinDepth = 5;
        constexpr uint32_t MaxDepth = 100;
        Color3f throughput(1.0f);
        Color3f L(0.0f);
        Ray3f ray(ray_);
        uint32_t depth=1;
        Intersection its;
        while(depth<MaxDepth && scene->rayIntersect(ray,its))
        {
            // std::cout<<"line24"<<std::endl;
            if(its.mesh->isEmitter())
            {
                LightQueryRecord record(its.p,ray.o,its.shFrame.n);
                L+=throughput*its.mesh->getEmitter()->eval(its.mesh,record);
            }

            //Russian Roullete
            if(depth>MinDepth)
            {
                float p=std::min(throughput.maxCoeff(),0.99f);
                if(sampler->next1D() > p)
                    break;
                throughput/=p;
            }
            
            //sample directlight
            auto meshes = scene->getMeshes();
            Color3f Ld(0.f);

            for(auto mesh:meshes)
            {
                if(mesh && mesh->isEmitter())
                {
                    if(mesh == its.mesh)
                        continue;
                    
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
                            Ld += radiance * its.mesh->getBSDF()->eval(bRec) * std::abs(its.shFrame.n.dot(wi)) * std::abs(record.AreaLight_normal.dot(-wi))
                            /(record.Light_Sample_point - its.p).squaredNorm()/pdf;
                        }
                    }
                    
                }
            }
            L+=throughput*Ld;



            //sample indirectlight
            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
            Color3f bsdf=its.mesh->getBSDF()->sample(bRec,sampler->next2D());
            throughput*=bsdf;
            depth++;
            ray=Ray3f(its.p,its.shFrame.toWorld(bRec.wo));

        }
    
        return L/2.f;
    #endif
    
    

    }

    std::string toString() const {
        return "PathEms[]";
    }
};

NORI_REGISTER_CLASS(PathEms, "path_ems");

NORI_NAMESPACE_END
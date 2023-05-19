#include<nori/integrator.h>
#include<nori/scene.h>
#include<nori/emitter.h>
#include<nori/sampler.h>
#include<nori/warp.h>
#include<nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMats : public Integrator {
    public:
    PathMats(const PropertyList &props) {}

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_)const{
    #if 1
        constexpr uint32_t MinDepth = 5;
        constexpr uint32_t MaxDepth = 1000;
        Color3f throughput(1.0f);
        Color3f L(0.0f);
        Ray3f ray(ray_);
        uint32_t depth=1;
        Intersection its;
        while(depth<MaxDepth && scene->rayIntersect(ray,its))
        {
            if(its.mesh->isEmitter())
            {
                EmitterQueryRecord record(ray.o,its.p,its.shFrame.n);
                L+=throughput*its.mesh->getEmitter()->eval(record);
            }

            //Russian Roullete
            if(depth>MinDepth)
            {
                float p=std::min(throughput.maxCoeff(),0.99f);
                if(sampler->next1D() > p)
                    break;
                throughput/=p;
            }
        

            //sample indirect light
            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
            Color3f bsdf=its.mesh->getBSDF()->sample(bRec,sampler->next2D());
            throughput*=bsdf;
            depth++;
            ray=Ray3f(its.p,its.shFrame.toWorld(bRec.wo));

        }
    
        return L;
    #endif
    
    

    }

    std::string toString() const {
        return "PathMats[]";
    }
};

NORI_REGISTER_CLASS(PathMats, "path_mats");

NORI_NAMESPACE_END
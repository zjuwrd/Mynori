#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {
 public:
  PathMisIntegrator(const PropertyList& props) {}

  Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
    constexpr uint32_t Max_depth = 1000;
    constexpr uint32_t Min_depth = 5;
    Color3f L(0.0f);
    Color3f throughput(1.f);
    Ray3f current_ray(ray);
    float probablity = 1.f;
    float w_mats = 1.f;
    int depth = 1;
    Intersection its;
    if(!scene->rayIntersect(current_ray,its))
    {
        return L;
    }
    
    while(depth<Max_depth)
    {
        if(its.mesh->isEmitter())
        {
            EmitterQueryRecord LRec(its.p,current_ray.o,its.shFrame.n);
            L += w_mats*throughput*its.mesh->getEmitter()->eval(LRec);
            if(std::isnan(L[0])||std::isnan(L[1])||std::isnan(L[2]))
            {
                L=0.f;
            }
        }

        float light_pdf = 0.f;
        const Mesh* light = scene->SampleLight(sampler->next1D(),light_pdf);
;        
        EmitterQueryRecord LRec(its.p);
        float em_pdf;
        //direct light
        Color3f Ld = light->getEmitter()->sample(LRec,sampler);
        em_pdf = LRec.pdf;
        auto wi = LRec.wi;
        Ray3f shadow_ray(its.p,wi,Epsilon,LRec.dist-Epsilon);  

        if(!scene->rayIntersect(shadow_ray))//no occlusion
        {
            float cosTheta = std::max(0.f, its.shFrame.n.dot(wi));
            BSDFQueryRecord BRec(its.shFrame.toLocal(-current_ray.d),its.shFrame.toLocal(wi),ESolidAngle);
            Color3f f = its.mesh->getBSDF()->eval(BRec);
            float mat_pdf = its.mesh->getBSDF()->pdf(BRec);
            float w_ems = mat_pdf+em_pdf>0.f? em_pdf/(mat_pdf+em_pdf):0.f;
            L += Ld*f*cosTheta*w_ems*throughput;
            if((L[0])<0.f||(L[1])<0.f||(L[2])<0.f)
            {
                std::cout<<"line67:L is nan"<<std::endl;
               
                L=0.f;
            }
        }

        //Russian Rullete
        if(depth>Min_depth)
        {
            float p = std::min(throughput.maxCoeff(),0.99f);
            if(sampler->next1D()>p)
            {
                break;
            }
            throughput /= p;
        }

        //indirect light
        BSDFQueryRecord BRec(its.shFrame.toLocal(-current_ray.d));
        Color3f f = its.mesh->getBSDF()->sample(BRec,sampler->next2D());
        throughput*=f;
        
        float mat_pdf = its.mesh->getBSDF()->pdf(BRec);
        current_ray = Ray3f(its.p,its.toWorld(BRec.wo));
        if(!scene->rayIntersect(current_ray,its))
        {
            break;
        }

        
        if(BRec.measure == EDiscrete)
        {
            w_mats = 1.f;
        }
        else if(its.mesh->isEmitter() )
        {
            EmitterQueryRecord LRec(current_ray.o,its.p,its.shFrame.n);
            float new_pdf_em = its.mesh->getEmitter()->pdf(LRec);
            w_mats = mat_pdf + new_pdf_em > 0.f? mat_pdf/(mat_pdf+new_pdf_em):0.f;
        
        }
        
        
        ++depth;


    }

    if(std::isnan(L[0])||std::isnan(L[1])||std::isnan(L[2]))
            {
                // std::cout<<"line116:L is nan"<<std::endl;
                L=0.f;
            }

    return L;
  
    static int count = 0;
    if(count%1000==0)
    {
        std::cout<<"path "<<count<<"in mis"<<std::endl;
    }


  }

  std::string toString() const {
    return "PathMisIntegrator[]";
  }
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END
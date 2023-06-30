#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathEms : public Integrator {
 public:
  PathEms(const PropertyList& props) {}


#if 1
  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;           
        Ray3f iterRay(ray);         
        Color3f throughput = 1.f;   
        Color3f fr;                 
        Color3f L = 0.f;             
        bool lastSpecular = false;

        for (uint32_t depth = 0;scene->rayIntersect(iterRay, its) ; ++depth) {
            
            if ((depth == 0 || lastSpecular) && its.mesh->isEmitter()) {
                EmitterQueryRecord eQ = EmitterQueryRecord(iterRay.o,its.p,its.shFrame.n);
                L += throughput * its.mesh->getEmitter()->eval(eQ);
            }
            
            //sample direct light
            if(!its.mesh->getBSDF()->isDelta())
            {
              L += throughput * scene->SampleLd(iterRay,its,sampler);
              lastSpecular = false;
            }
            else
            {
              lastSpecular = true;
            }

            //Russian roulette
            float p = std::min(0.99f, throughput.maxCoeff());
            if (depth > 3) {
                if (sampler->next1D() > p) {
                    break;
                }
                throughput /= p;
            }

            //sample direction for indirect light
            BSDFQueryRecord bQ(its.shFrame.toLocal(-iterRay.d).normalized());
            fr = its.mesh->getBSDF()->sample(bQ, sampler->next2D());
            throughput *= fr;
            if (fr.getLuminance() == 0.0f) break;
            //generate new ray
            iterRay = Ray3f(its.p, its.toWorld(bQ.wo));
            
        }
        return L;
    }



#endif
  std::string toString() const {
    return "PathEms[]";
  }
};

NORI_REGISTER_CLASS(PathEms, "path_ems");
NORI_NAMESPACE_END
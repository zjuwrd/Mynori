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
  Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& _ray) const override {
#if 0
    
    Color3f L(0.f);
    Intersection its;
    auto w_ems = 0.5f;
    auto w_mat = 0.5f;
    
    
    if(!scene->rayIntersect(_ray,its))
    {
      return L;
    }

    auto mesh = its.mesh;
    //Le
    if(mesh->isEmitter())
    {
      EmitterQueryRecord eQE(_ray.o, its.p, its.shFrame.n);
      L += mesh->getEmitter()->eval(eQE);
    }

    //ld
    auto Rand_Emitter = scene->SampleLight(sampler->next1D())->getEmitter();
    EmitterQueryRecord eQ(its.p);
    Color3f Ld = Rand_Emitter->sample(eQ, sampler);
    Color3f bsdf = mesh->getBSDF()->eval(BSDFQueryRecord( its.toLocal( -_ray.d) , its.toLocal(eQ.wi), EMeasure::ESolidAngle));
    float Costheta = std::max(0.f, Frame::cosTheta(its.toLocal(eQ.wi)));
    
    if(!scene->rayIntersect(eQ.shadowRay())) // no occlution
      L += w_ems * Ld * bsdf * Costheta *scene->getEmissiveMeshesCount();

    //russian_roulette
    float RR_prob= std::max(0.001f, std::min(0.99f, _ray.throughput));
    if(sampler->next1D()>RR_prob)
    {
      return L;
    }



    //indirect
    BSDFQueryRecord bRec(its.toLocal(-_ray.d));
    Color3f bsdf_cosTheta = mesh->getBSDF()->sample(bRec, sampler->next2D());
    Ray3f ray_indirect(its.p, its.toWorld(bRec.wo));
    ray_indirect.depth = _ray.depth+1;
    ray_indirect.throughput *= (bsdf_cosTheta.maxCoeff() / RR_prob);

    Color3f L_indirect = Li(scene, sampler, ray_indirect);
    w_mat = its.mesh->getBSDF()->isDelta() ? 1.f:0.5f;

    L += w_mat * L_indirect * bsdf_cosTheta / RR_prob;

  




  return L;
#endif


#if 1
    constexpr size_t Max_trace_depth = 1000;
    constexpr size_t Min_trace_depth = 5;
    
    Color3f color = 0;
    Color3f throughput = 1;
    Ray3f rayRecursive = _ray;
    float probability;
    size_t depth = 0;
    Intersection its;
    bool is_specular=false;
    int cnt = 0;
    cnt++;
    
    while (scene->rayIntersect(rayRecursive, its) && depth<Max_trace_depth) {
      //Le
      if (its.mesh->isEmitter() && (is_specular || depth == 0)) {
        EmitterQueryRecord eQE(rayRecursive.o, its.p, its.shFrame.n);
        color += throughput * its.mesh->getEmitter()->eval(eQE);
        if(!color.isValid())
        {
          color = 0.f;
        }
      }

      //Ld
      if(!its.mesh->getBSDF()->isDelta())
      {  
        float light_pdf=0.f;
        const Emitter* light = scene->SampleLight(sampler->next1D(),light_pdf)->getEmitter();
        EmitterQueryRecord eQ(its.p);
        Color3f Ld = light->sample(eQ, sampler);
        Color3f bsdf = its.mesh->getBSDF()->eval(BSDFQueryRecord( its.toLocal(eQ.wi), its.toLocal(-rayRecursive.d), EMeasure::ESolidAngle));
        float cosTheta = Frame::cosTheta(its.toLocal(eQ.wi));
        if(!scene->rayIntersect(eQ.shadowRay())) // no occlusion
          color += throughput * bsdf * std::max(0.f,  cosTheta) * Ld / light_pdf;
        is_specular = false;
      }
      else{
        is_specular = true;
      }

      //russian roulete
      float p = std::min(throughput.maxCoeff(),0.99f);
      if(sampler->next1D()>p)
      {
        break;
      }
      else{
        throughput/=p;
      }
      
      //indirect light
      BSDFQueryRecord brec(its.toLocal(-rayRecursive.d));
      Color3f bsdf_cosTheta = its.mesh->getBSDF()->sample(brec,sampler->next2D());
      throughput *= bsdf_cosTheta;
      rayRecursive = Ray3f(its.p,its.toWorld(brec.wo));
      depth++;

    }

    return color;

#endif
  }
#endif

#if 0
  Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;           // Intersection information
        Ray3f iterRay(ray);         // Ray used for path tracing
        Color3f throughput = {1};   // Path throughput
        Color3f fr;                 // Reflection or refraction factor (albedo)
        Color3f L = {0};            // Final rendered light 
        std::vector<Mesh *> meshes; // All meshes of the scene
        EmitterQueryRecord eQ(Point3f(0.f,0.f,0.f));      // The emitter query record
        float eta = 1.0f;           // Cumulative eta
        float prob;                 // Probability to continue
        bool foundIntersection;     // Whether there is an intersection point found
        bool specularBounce = false;// Whether last bounce is a specular bounce

        /* Start the path tracing */
        for (uint32_t bounces = 0;scene->rayIntersect(iterRay, its) ; ++bounces) {
            /* Add the emitter light at path vertex */
            if ((bounces == 0 || specularBounce) && its.mesh->isEmitter()) {
                eQ = EmitterQueryRecord(iterRay.o,its.p,its.shFrame.n);
                L += throughput * its.mesh->getEmitter()->eval(eQ);
            }

            // if (its.mesh->getBSDF()->isNull()) {
            //     Ray3f temp = Ray3f(its.p, iterRay.d);
            //     memcpy(&iterRay, &temp, sizeof(Ray3f));
                
            //     bounces--;
            //     continue;
            // }
            
            /* Explicit sampling direct emitters */
            L += throughput * scene->SampleLd(iterRay,its,sampler);
            
            
            /* Account for indirect light, we sample a new direction on this surface */
            BSDFQueryRecord bQ(its.shFrame.toLocal(-iterRay.d).normalized());
            fr = its.mesh->getBSDF()->sample(bQ, sampler->next2D());
            specularBounce = bQ.measure == EDiscrete;
            if (fr.getLuminance() == 0.0f) break;

            /* Update the iteration ray given the sampling */
            // Ray3f temp = Ray3f(its.p, its.shFrame.toWorld(bQ.wo));
            // memcpy(&iterRay, &temp, sizeof(Ray3f));
            iterRay.o = its.p;
            iterRay.d = its.toWorld(bQ.wo);
            iterRay.update();
            // iterRay = Ray3f(its.p, its.toWorld(bQ.wo));


            /* Compute the probability to continue */
            eta *= bQ.eta;
            prob = std::min(0.99f, throughput.maxCoeff() * eta * eta);
            /* Only start doing Russian Roulette after at least three bounces */
            prob = bounces < 4 ? 1 : prob;

            /* Use the Russian Roulette */
            if (sampler->next1D() < prob) {
                throughput *= fr / prob;
            } else {
                break;
            }
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
#if 1

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {
 private:
 // Estimate Direct Lighting using MIS
	// return appropriately weighted terms.
	Color3f LiDirect(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& isect) const
	{
		Color3f L_ems(0.0f), L_mats(0.0f);		
		const BSDF* bsdf = isect.mesh->getBSDF();
		
		Color3f dbg(0.f);

		// Choose a light
		float pdf = 1.0f / scene->getLights().size();
		const Emitter* random_emitter = scene->SampleLight(sampler->next1D())->getEmitter();

		// Emitter Sampling
		// Perform only if not a delta bsdf		
		if(!bsdf->isDelta())
		{
			EmitterQueryRecord eRec(isect.p);
			float pdf_e, pdf_m;
			

			Color3f Li = random_emitter->sample(eRec, sampler);
			pdf_e = eRec.pdf;
			
			BSDFQueryRecord bRec(isect.toLocal(-ray.d), isect.toLocal(eRec.wi), ESolidAngle);
			Color3f f = bsdf->eval(bRec);
			pdf_m = bsdf->pdf(bRec);
			if (pdf_e != 0.0f)
			{
				float mis = pdf_e / (pdf_m + pdf_e);

				// Compute lighting
				L_ems = f * Li * fabsf(isect.shFrame.n.dot(eRec.wi));

				// Compute shadow ray only when 
				if (L_ems.isValid() && !L_ems.isZero())
				{
					// Trace a shadow ray only now
					float V = scene->rayIntersect(Ray3f(isect.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist)) ? 0.0f : 1.0f;
					L_ems *= V;
				}
				else
					L_ems = Color3f(0.0f);

				dbg = L_ems;

				if (!random_emitter->isDelta())
				{
					// The BSDF has no way of generating a direction that would hit this light
					// Hence multiply by MIS value only when 
					L_ems *= mis;
				}				
			}			
		}
				
		// BSDF sampling
		// If the chosen light was a delta light, we can't expect the BSDF to sample a direction that will hit the light
		dbg = 0.f;
		if (!random_emitter->isDelta())
		{
			BSDFQueryRecord bRec(isect.toLocal(-ray.d));
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			float pdf_m = bsdf->pdf(bRec);

			if (!f.isZero() && pdf_m != 0.0f && !isnan(pdf_m))
			{
				Intersection light_isect;
				if (scene->rayIntersect(Ray3f(isect.p, isect.toWorld(bRec.wo), Epsilon, INFINITY), light_isect))
				{
					// check if a light soruce
					if (light_isect.mesh->isEmitter() && light_isect.mesh->getEmitter() == random_emitter)
					{
						const Emitter* light = light_isect.mesh->getEmitter();

						EmitterQueryRecord eRec(isect.p,light_isect.p,light_isect.shFrame.n);
						
						Color3f Li = light->eval(eRec);
						float mis = 1.0f;

						// If the BSDF was delta, the light could have never generated this reverse direction
						if (!bsdf->isDelta())
						{
							float pdf_e = light->pdf(eRec);
							mis = pdf_m / (pdf_m + pdf_e);
						}

						L_mats = (f * Li * fabsf(Frame::cosTheta(bRec.wo)));
						dbg = L_mats;
						L_mats *= mis;
					}
				}
			}
		}
		
		// return dbg/pdf;
		// Divide by the pdf of choosing the random light
		return (L_ems + L_mats) / pdf;
	}

 
 public:
  PathMisIntegrator(const PropertyList& props) {}

  Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_) const override {

#if 0
	 Color3f L(0.f);
    Intersection its;
    
    if(!scene->rayIntersect(ray_,its))
    {
      return L;
    }

    auto mesh = its.mesh;
    //Le
    if(mesh->isEmitter()&&ray_.depth == 0)
    {
      EmitterQueryRecord lRecE(ray_.o, its.p, its.shFrame.n);
      L += mesh->getEmitter()->eval(lRecE);
    }

    //ld
	auto Ld = this->LiDirect(scene,sampler,ray_,its);
	L += Ld/2.f;


    //russian_roulette
	
    float RR_prob=1.f;
	if(ray_.depth > 4)
	{
		RR_prob = std::max(0.001f, std::min(0.99f, ray_.throughput));
		if(sampler->next1D()>RR_prob)
		{
			return L;
		}
	}
    //indirect
    BSDFQueryRecord bRec(its.toLocal(-ray_.d));
    Color3f bsdf_cosTheta = mesh->getBSDF()->sample(bRec, sampler->next2D());
    Ray3f ray_indirect(its.p, its.toWorld(bRec.wo));
    ray_indirect.depth = ray_.depth+1;
    ray_indirect.throughput *= (bsdf_cosTheta.maxCoeff() / RR_prob);

    Color3f L_indirect = Li(scene, sampler, ray_indirect);
    if(!its.mesh->getBSDF()->isDelta())
	{
		L_indirect *= 0.5f;
	}

    L += L_indirect * bsdf_cosTheta / RR_prob;

  return L;


#endif



#if 1
    constexpr size_t Max_depth = 1000;
    constexpr size_t Min_depth = 5;

    Color3f throughput = 1.f;
    Color3f L = 0.f;
    Ray3f curRay(ray_);
    Intersection its;
    size_t depth = 0;

    while(scene->rayIntersect(curRay, its)&& depth<Max_depth)
    {
      //Le
      if(its.mesh->isEmitter() && depth == 0)
      {
        EmitterQueryRecord ERec(curRay.o,its.p,its.shFrame.n);
        L += throughput*its.mesh->getEmitter()->eval(ERec);
      }
      //sample direct light using MIS
      Color3f Ld=0.f;
      Ld = this->LiDirect(scene,sampler,curRay,its);

      L += Ld*throughput /2.f;

      //Russian Roulete
      if(depth > Min_depth)
      {
        float p = std::min(throughput.maxCoeff(),0.99f); 
        if(sampler->next1D() > p)
        {
          break;
        }
        throughput/=p;
      }
      
	  //sample a new ray
      {
        BSDFQueryRecord BRec(its.toLocal(-curRay.d));
        Color3f bsdf = its.mesh->getBSDF()->sample(BRec,sampler->next2D());
        throughput *= bsdf;
        curRay = Ray3f(its.p,its.toWorld(BRec.wo));
		if(!its.mesh->getBSDF()->isDelta())
		{
			throughput *= 0.5f;
		}

	  }

      depth++;

    }
    return L;
#endif


  }

  std::string toString() const {
    return "PathMisIntegrator[]";
  }
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END

#endif
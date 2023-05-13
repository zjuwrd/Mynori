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

#if 1
	constexpr uint32_t MAX_depth = 40;
	constexpr uint32_t Min_depth = 3; 
	Color3f Li(0.f);
	Color3f throughput(1.f);
	Intersection its;
	Ray3f ray(ray_);
	float w_mats=1.f;
	float w_ems = 1.f;



	for(uint32_t depth=0;depth<MAX_depth && scene->rayIntersect(ray,its); ++depth)
	{
		//Le
		if(its.mesh->isEmitter())
		{
			EmitterQueryRecord eQ(ray.o,its.p,its.shFrame.n);
			const Emitter* emitter = its.mesh->getEmitter();
			Li+= w_mats * throughput * emitter->eval(eQ);
		}

		//Ld
		if(its.mesh->getBSDF()->isDelta())
		{
			//do nothing
		}
		else{
			EmitterQueryRecord eQ(its.p);
			float light_pdf =0.f;
			const Emitter* emitter = scene->SampleLight(sampler->next1D(),light_pdf)->getEmitter();
			Color3f Ld = emitter->sample(eQ,sampler);
			float ems_pdf = eQ.pdf;
			BSDFQueryRecord bQ(its.toLocal(-ray.d),its.toLocal(eQ.wi),ESolidAngle);
			Color3f fr = its.mesh->getBSDF()->eval(bQ);
			float mats_pdf = its.mesh->getBSDF()->pdf(bQ);
			float cosTheta = std::max(0.f, Frame::cosTheta(its.toLocal(eQ.wi)));
			w_ems = ems_pdf+mats_pdf>0.f?ems_pdf/(ems_pdf+mats_pdf):0.f;

			Li += w_ems * throughput * fr * Ld * cosTheta / light_pdf;
		}
	

		//russian roulette
		if(depth>Min_depth)
		{
			float q = std::max(0.05f,1.f-throughput.maxCoeff());
			if(sampler->next1D()<q)
				break;
			throughput /= 1.f-q;
		}

		//sample indirect light
		BSDFQueryRecord bQ(its.toLocal(-ray.d));
		Color3f bsdf_cosTheta = its.mesh->getBSDF()->sample(bQ,sampler->next2D());
		throughput *= bsdf_cosTheta;

		if(bQ.measure == EDiscrete)
		{
			w_mats = 1.f;
		}
		else{
			float mats_pdf = its.mesh->getBSDF()->pdf(bQ);
			Ray3f testray(its.p,its.toWorld(bQ.wo));
			Intersection testits;
			if(scene->rayIntersect(testray,testits))
			{
				if(testits.mesh->isEmitter())
				{
					EmitterQueryRecord eQ(its.p,testits.p,testits.shFrame.n);
					const Emitter* emitter = testits.mesh->getEmitter();
					float ems_pdf = emitter->pdf(eQ);
					w_mats = mats_pdf+ems_pdf>0.f?mats_pdf/(mats_pdf+ems_pdf):0.f;
				}
				else
				{
					w_mats = 1.f;
				}
			}
			else
			{
				w_mats = 1.f;
			}
		}

		ray = Ray3f(its.p,its.toWorld(bQ.wo));
	}

	

	return Li;
	 


#endif



#if 0
    constexpr size_t Max_depth = 1000;
    constexpr size_t Min_depth = 5;

    Color3f throughput = 1.f;
    Color3f L = 0.f;
    Ray3f curRay(ray_);
    Intersection its;
    size_t depth = 0;

	bool is_specular = false;

    while(scene->rayIntersect(curRay, its)&& depth<Max_depth)
    {
      //Le
      if(its.mesh->isEmitter() && (depth == 0 || is_specular))
      {
        EmitterQueryRecord ERec(curRay.o,its.p,its.shFrame.n);
        L += throughput*its.mesh->getEmitter()->eval(ERec);
      }
      //sample direct light using MIS
      Color3f Ld=0.f;
	  Ld = this->LiDirect(scene,sampler,curRay,its);

	
      L += Ld*throughput;

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
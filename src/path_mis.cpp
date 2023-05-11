#if 0
#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include<nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PathIntegratorMis : public Integrator
{
public:

	enum class DirectSamplingStrategy
	{
		SAMPLE_ALL_LIGHTS,
		SAMPLE_ONE_LIGHT
	};

	PathIntegratorMis(const PropertyList& props)
	{
		m_rrStart = props.getInteger("rrStart", 5);
		m_maxDepth = props.getInteger("maxDepth", -1);
		
		// Which strategy to use for our NEE scheme
		std::string strategy = props.getString("strategy", "sample_one_light");
		if (strategy == "sample_all_lights")
			m_strategy = DirectSamplingStrategy::SAMPLE_ALL_LIGHTS;
		else m_strategy = DirectSamplingStrategy::SAMPLE_ONE_LIGHT;

	}

	~PathIntegratorMis() {}

	// Estimate Direct Lighting using MIS
	// return appropriately weighted terms.
	Color3f LiDirect(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& isect) const
	{
		Color3f L_ems(0.0f), L_mats(0.0f);		
		const BSDF* bsdf = isect.mesh->getBSDF();

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
						L_mats *= mis;
					}
				}
			}
		}
		
		// Divide by the pdf of choosing the random light
		return (L_ems + L_mats) / pdf;
	}


	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		Intersection isect;

		Color3f L(0.0f);
		Color3f throughput(1.0f);
		int depth = 0;
		Ray3f traced_ray = ray;
		bool wasLastBounceSpecular = false;

		while (depth < m_maxDepth || m_maxDepth == -1)
		{
			// Check if ray misses the scene
			if (!scene->rayIntersect(traced_ray, isect))
			{
				L += throughput *Color3f(0.f);
				break;
			}

			// Check if direct emitter intersection
			// return this only if direct hit or previous hit was from a specular surface
			// because we couldnt have sampled it using a light sampling strategy.
			if (isect.mesh->isEmitter() && (depth == 0))
			{
				//if (depth == 0 || wasLastBounceSpecular)
				{
					EmitterQueryRecord eRec(traced_ray.o,isect.p,isect.shFrame.n);
					eRec.wi = traced_ray.d;
					eRec.n = isect.shFrame.n;
					L += throughput * isect.mesh->getEmitter()->eval(eRec);
				}
				//else break;

				// Assume for now we dont bounce off the light sources.				
			}
			

			const BSDF* bsdf = isect.mesh->getBSDF();			

			// NEE
			Color3f Li = LiDirect(scene, sampler, traced_ray, isect);
			Color3f debug = throughput * Li;
			L += throughput * Li;
						
			// Sample a reflection ray
			BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			Vector3f reflected_dir = isect.toWorld(bRec.wo);

			throughput *= f * fabsf(Frame::cosTheta(bRec.wo));

			// Check if specular bounce
			wasLastBounceSpecular = bsdf->isDelta();

			// Check if we've reached a zero throughput. No point in proceeding further.
			if (throughput.isZero())
				break;

			// Check for russian roulette
			if (depth > m_rrStart)
			{
				float rrprob = throughput.getLuminance();
				if (sampler->next1D() > rrprob)
					break;
				else throughput /= rrprob;
			}
			else if (depth > m_maxDepth && m_maxDepth != -1)
			{
				// forcibly terminate
				break;
			}

			// Propogate
			traced_ray = Ray3f(isect.p, reflected_dir, Epsilon, INFINITY);
			depth++;
		}
		
		return L;
	}

	std::string toString() const
	{
		return tfm::format("PathIntegratorMis[\nrrStart = %d\n]", m_rrStart);
	}

private:
	int m_rrStart;				// from which bounce should russian roulette start.
	int m_maxDepth;				// Fixed length cutoff
	DirectSamplingStrategy m_strategy;

};

NORI_REGISTER_CLASS(PathIntegratorMis, "path_mis")
NORI_NAMESPACE_END

#endif

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
						L_mats *= mis;
					}
				}
			}
		}
		
		// Divide by the pdf of choosing the random light
		return (L_ems + L_mats) / pdf;
	}

 
 public:
  PathMisIntegrator(const PropertyList& props) {}

  Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_) const override {
    constexpr size_t Max_depth = 1000;
    constexpr size_t Min_depth = 5;

    #if 0
    
    Color3f L = 0.f;
    Color3f throughput = 1.f;
    Ray3f curRay = ray_;
    float RR_prob=0.f;
    float w_mats = 1.0f,w_ems = 1.0f;
    int depth = 0;
    Intersection its;
    if (!scene->rayIntersect(curRay, its)) {
      return L;
    }
    while (true) {

      //Le
      if(its.mesh->isEmitter())
      {
        EmitterQueryRecord ERec(curRay.o,its.p,its.shFrame.n);
        L += w_mats*throughput*its.mesh->getEmitter()->eval(ERec);
      }

      //Ld
      {
        float LightPdf=0.f;
        auto LightMesh = scene->SampleLight(sampler->next1D(),LightPdf);
        EmitterQueryRecord ERec(its.p);
        Color3f Ld = LightMesh->getEmitter()->sample(ERec,sampler);
        if(!scene->rayIntersect(ERec.shadowRay()))
        {  

          BSDFQueryRecord BRec(its.toLocal(ERec.wi),its.toLocal(-curRay.d),EMeasure::ESolidAngle);
          Color3f bsdf = its.mesh->getBSDF()->eval(BRec);
          float mats_pdf = its.mesh->getBSDF()->pdf(BRec);
          float ems_pdf = LightPdf * ERec.pdf;
          w_ems = (mats_pdf+ems_pdf>0.f)?ems_pdf/(mats_pdf+ems_pdf):0.f;
          float abscosTheta = std::abs(Frame::cosTheta(its.toLocal(ERec.wi)));
          L += w_ems*throughput*bsdf*abscosTheta* Ld;
      
        }
      
      }

      //Russian Roulete
      if(depth>Min_depth)
      {
        RR_prob = std::min(throughput.maxCoeff(),0.99f);
        if(sampler->next1D() > RR_prob)
        {
          break;
        }
        else{
          throughput/=RR_prob;
        }
      }      

      //generate new Ray
      {
        BSDFQueryRecord BRec(its.toLocal(-curRay.d));
        Color3f bsdf = its.mesh->getBSDF()->sample(BRec,sampler->next2D());
        throughput *= bsdf;
        curRay = Ray3f(its.p,its.toWorld(BRec.wo));
        
        if(!scene->rayIntersect(curRay,its)){break;}   
      
        //generating new w_mats
        if(its.mesh->isEmitter())
        {
            if(BRec.measure == EDiscrete)
            {
              w_mats = 1.f;
            }
            else
            {
              EmitterQueryRecord eRec(curRay.o,its.p,its.shFrame.n);
              float ems_pdf = its.mesh->getEmitter()->pdf(eRec);
              float mats_pdf = its.mesh->getBSDF()->pdf(BRec);
              w_mats = (mats_pdf+ems_pdf>0.f)?mats_pdf/(mats_pdf+ems_pdf):0.f;
            }
        }
      }        
      
      depth++;
    }

    return L;
    #endif
    
    #if 1
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
      // {
      //   float w_ems=0.f,w_mat=0.f;
      //   Color3f L_ems=0.f,L_mat=0.f;
        
      //   auto EmisMesh = scene->SampleLight(sampler->next1D());
      //   const auto emitter = EmisMesh->getEmitter();
            
      //   //sample from light
      //   {
      //     //only do this if material is not perfect specular, 
      //     //in this case, the material can only be diffuse
      //     if(its.mesh->getBSDF()->isDiffuse())
      //     {  
      //       EmitterQueryRecord ERec(its.p);
      //       Color3f Ld = emitter->sample(ERec,sampler);
      //       //check if the light is visible
      //       if(!scene->rayIntersect(ERec.shadowRay()))
      //       {
      //         BSDFQueryRecord BRec(its.toLocal(ERec.wi),its.toLocal(-curRay.d),EMeasure::ESolidAngle);
      //         Color3f bsdf = its.mesh->getBSDF()->eval(BRec);
              
      //         float mats_pdf = its.mesh->getBSDF()->pdf(BRec);
      //         float ems_pdf = emitter->pdf(ERec);
              
      //         w_ems = (mats_pdf+ems_pdf>0.f)?ems_pdf/(mats_pdf+ems_pdf):0.f;
      //         float abscosTheta = std::abs(Frame::cosTheta(its.toLocal(ERec.wi)));
      //         L_ems = bsdf*abscosTheta* Ld;
      //       }
      //     }
      //   }
      //   //std::cout<<"line 277"<<std::endl;
      //   //sample from BSDF
      //   {
      //     BSDFQueryRecord BRec(its.toLocal(-curRay.d));  
      //     Color3f bsdf_abscosTheta = its.mesh->getBSDF()->sample(BRec,sampler->next2D());
          
      //     Ray3f testRay(its.p,its.toWorld(BRec.wo));
      //     Intersection testIts;
      //     //std::cout<<"line 285"<<std::endl;
      //     if(scene->rayIntersect(testRay,testIts) )
      //     {
      //       if(testIts.mesh == nullptr)
      //       {
      //         //std::cout<<"testIts.esh == nullptr"<<std::endl;
      //       }
      //       if(testIts.mesh->getEmitter() == emitter)
      //       {  
      //         //std::cout<<"line 288"<<std::endl;
      //         EmitterQueryRecord ERec(its.p,testIts.p,testIts.shFrame.n);
      //         Color3f Ld = emitter->eval(ERec);
      //         //std::cout<<"line 291"<<std::endl;
      //         float mats_pdf = its.mesh->getBSDF()->pdf(BRec);
      //         float ems_pdf = emitter->pdf(ERec);
      //         //std::cout<<"line 294"<<std::endl;
      //         if(its.mesh->getBSDF()->isDiffuse())
      //           w_mat = (mats_pdf+ems_pdf>0.f)?mats_pdf/(mats_pdf+ems_pdf):0.f;
      //         else  //specular sitution
      //           w_mat = 1.f; 
      //         //std::cout<<"line 298"<<std::endl;

      //         L_mat = bsdf_abscosTheta*Ld;   
      //       }
      //     }
      //   }
      //   Ld = L_ems*w_ems + L_mat*w_mat;
      // }
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
      //std::cout<<"line 318"<<std::endl;
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
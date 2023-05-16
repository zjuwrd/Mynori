#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator { 
 public:
  PathMisIntegrator(const PropertyList& props) {}

  Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_) const override {

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
			Li+= throughput * emitter->eval(eQ);
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
			if(!scene->rayIntersect(eQ.shadowRay())) // no occlusion
			{
				float ems_pdf = eQ.pdf;
				BSDFQueryRecord bQ(its.toLocal(-ray.d),its.toLocal(eQ.wi),ESolidAngle);
				Color3f fr = its.mesh->getBSDF()->eval(bQ);
				float mats_pdf = its.mesh->getBSDF()->pdf(bQ);
				float cosTheta = std::max(0.f, Frame::cosTheta(its.toLocal(eQ.wi)));
				w_ems = ems_pdf+mats_pdf>0.f?ems_pdf/(ems_pdf+mats_pdf):0.f;

				Li += w_ems * throughput * fr * Ld * cosTheta / light_pdf;
			}
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
		Color3f bsdf_cosTheta = its.mesh->getBSDF()->sample(bQ, sampler);
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
		throughput *= w_mats;

		ray = Ray3f(its.p,its.toWorld(bQ.wo));
	}

	return Li;
  }

  std::string toString() const {
    return "PathMisIntegrator[]";
  }
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END


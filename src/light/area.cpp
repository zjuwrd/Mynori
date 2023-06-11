#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter {
public:
    AreaEmitter(const PropertyList &props) :m_radiance(props.getColor("radiance")){
		m_type = EmitterType::EMITTER_AREA;
    }

    virtual std::string toString() const {
        return tfm::format(
                "AreaLight[\n"
                "  radiance = %s,\n"
                "]",
                m_radiance.toString());
    }

    virtual Color3f getRadiance() const {return m_radiance;}

	virtual Color3f eval(const EmitterQueryRecord & lRec) const {
        assert(m_mesh);
        
		if (lRec.n.dot(lRec.wi) < 0.0f)
			return m_radiance;
		else return 0.0f;
    }


    virtual Color3f sample(const Mesh* mesh, EmitterQueryRecord &eRec, Sampler* &sample) const override
    {
        MeshQeuryRecord point = mesh->UniformSamplePoint(sample);
        eRec.p = point.p;
        eRec.n = point.n;
        eRec.wi = (eRec.p - eRec.ref).normalized();
        eRec.pdf = pdf(mesh, eRec);

        if(eRec.pdf > 0.0f && !std::isnan(eRec.pdf) && !std::isinf(eRec.pdf)) {
            return eval(eRec) / eRec.pdf;
        }

        return Color3f(0.f);
    }

    virtual float pdf(const Mesh* mesh, const EmitterQueryRecord &eRec) const override
    {
        float pdf__=0.f;
        if(eRec.n.dot(-eRec.wi) >0.f) {
            pdf__ = mesh->getDpdf()->getNormalization(); 
        } 


        return pdf__;
    }


    virtual Color3f sample(EmitterQueryRecord & lRec, Sampler* sampler) const {
        //m_mesh shouldn't be null
        assert(m_mesh);
        
		// m_mesh->samplePosition(sampler->next2D(), lRec.p, lRec.n, sampler->next1D());
        MeshQeuryRecord res = m_mesh->UniformSamplePoint(sampler);
        lRec.p = res.p;
        lRec.n = res.n;
        lRec.wi = (lRec.p - lRec.ref).normalized();
        lRec.pdf = pdf(lRec);
        lRec.dist = std::sqrt( (lRec.p - lRec.ref).squaredNorm() );
		// fill 

		if(!std::isnan(lRec.pdf) && lRec.pdf > 0.0f && std::abs(lRec.pdf) != std::numeric_limits<float>::infinity()) 
        {
            return eval(lRec) / lRec.pdf;
        }
		else return 0.0f;
    }

    //for consistency, we use solid angle measure and pre-multiply the cosThetaPrime term
	virtual float pdf(const EmitterQueryRecord &lRec) const {
        assert(m_mesh);
        float pdf_ = 0.0f;

		float cosThetaPrime = std::abs(lRec.n.dot(-lRec.wi));

        //tranform pdf from area measure to solid angle measure    
        if(std::abs(cosThetaPrime)<1e-5f )
        {
            return 0.0f;
        }

		pdf_ = m_mesh->getDpdf()->getNormalization() * (lRec.p - lRec.ref).squaredNorm() / cosThetaPrime;
		//erase irrational values
        if (isnan(pdf_) || fabsf(pdf_) == INFINITY)
			pdf_ = 0.f;
	
    
    
    	return pdf_;
    }

    void setParent(NoriObject *parent)
	{
		auto type = parent->getClassType();
		if (type == EMesh)
			m_mesh = static_cast<Mesh*>(parent);
	}

    virtual PhotonRay ShootPhoton(Sampler* sampler)const{
        MeshQeuryRecord res = m_mesh->UniformSamplePoint(sampler);
        Vector3f dir_local = Warp::squareToCosineHemisphere(sampler->next2D());
        Frame frame(res.n);
        Vector3f dir_world = frame.toWorld(dir_local);

        Ray3f ray = Ray3f(res.p,dir_world,Epsilon,INFINITY);
        Color3f flux = M_PI * m_mesh->totalArea() * m_radiance;

        return PhotonRay(ray,flux);

    }

	

protected:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END
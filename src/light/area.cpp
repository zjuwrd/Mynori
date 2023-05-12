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

    virtual Color3f sample(EmitterQueryRecord & lRec, Sampler* sampler) const {
        //m_mesh shouldn't be null
        assert(m_mesh);
        
		m_mesh->samplePosition(sampler->next2D(), lRec.p, lRec.n, sampler->next1D());

		// fill in messages.
		lRec.wi = (lRec.p - lRec.ref).normalized();
		lRec.pdf = pdf(lRec);
        lRec.emitter = this;
		lRec.dist = (lRec.p - lRec.ref).norm();
		

		if(!std::isnan(lRec.pdf) && lRec.pdf > 0.0f && std::abs(lRec.pdf) != std::numeric_limits<float>::infinity()) 
        {
            return eval(lRec) / lRec.pdf;
        }
		else return 0.0f;
    }

    //for consistency, we use solid angle measure and pre-multiply the cosThetaPrime term
	virtual float pdf(const EmitterQueryRecord &lRec) const {
        assert(m_mesh);

		float cosThetaPrime = std::max(0.f , lRec.n.dot(-lRec.wi));

        //tranform pdf from area measure to solid angle measure    
        if(std::abs(cosThetaPrime)<1e-5f )
        {
            return 0.0f;
        }
		float pdf_ = m_mesh->pdf() * std::pow(lRec.dist,2) / cosThetaPrime;
		//erase irrational values
        if (isnan(pdf_) || fabsf(pdf_) == INFINITY)
			return 0.0f;
		return pdf_;
    }

    void setParent(NoriObject *parent)
	{
		auto type = parent->getClassType();
		if (type == EMesh)
			m_mesh = static_cast<Mesh*>(parent);
	}


	

protected:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END
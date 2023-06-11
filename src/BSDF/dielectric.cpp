/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include<nori/sampler.h>
NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        m_transmittance = propList.getFloat("transmittance", 1.f);

        m_eta = m_intIOR / m_extIOR;
        m_invEta = m_extIOR / m_intIOR; 
    }

    virtual bool isDelta()const{return true; }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, Sampler* sampler) const { return sample(bRec,sampler->next2D()); }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        
        #if 1
        bRec.measure = EDiscrete;
        float cosThetaI = Frame::cosTheta(bRec.wi);
        float kr = fresnel(cosThetaI, m_extIOR, m_intIOR);

        kr = FrDielectric(cosThetaI, m_eta);

        if (sample.x() < kr) {//reflect
        bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
        bRec.eta = 1.f;
        return Color3f(1.0f);
        } 
        else {//refract
            Vector3f n = Vector3f(0.0f, 0.0f, 1.0f);
            float factor = m_intIOR / m_extIOR;
            if (Frame::cosTheta(bRec.wi) < 0.f) {//transmit from outside to inside
                factor = m_intIOR / m_extIOR;
                n = -n;
                bRec.IOR_i = m_intIOR;
                bRec.IOR_o = m_extIOR;
            }
            else{//transmit from inside to outside
                factor = m_extIOR / m_intIOR;
                bRec.IOR_i = m_extIOR;
                bRec.IOR_o = m_intIOR;
            }
            
            if(!Refract(bRec.wi,n,m_eta,&bRec.eta,&bRec.wo)){
                bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
                bRec.eta = 1.f;
                return Color3f(1.0f);
            }
            
            // bRec.wo = refract(bRec.wi, n, factor);
            // bRec.eta = m_intIOR / m_extIOR;
            
            return Color3f(m_transmittance);
        }
        #endif



    }


    static inline float FrDielectric(float cosTheta_i, float eta) {
    cosTheta_i = std::clamp(cosTheta_i, -1.f, 1.f);
    // Potentially flip interface orientation for Fresnel equations
    if (cosTheta_i < 0) {
        eta = 1 / eta;
        cosTheta_i = -cosTheta_i;
    }

    // Compute $\cos\,\theta_\roman{t}$ for Fresnel equations using Snell's law
    float sin2Theta_i = 1 - (cosTheta_i*cosTheta_i);
    float sin2Theta_t = sin2Theta_i / (eta*eta);
    if (sin2Theta_t >= 1)
        return 1.f;
    float cosTheta_t = std::sqrt(1.f - sin2Theta_t);

    float r_parl = (eta * cosTheta_i - cosTheta_t) / (eta * cosTheta_i + cosTheta_t);
    float r_perp = (cosTheta_i - eta * cosTheta_t) / (cosTheta_i + eta * cosTheta_t);
    return ((r_parl*r_parl) + (r_perp * r_perp)) / 2;
}


    static bool Refract(Vector3f wi, Normal3f n, float eta, float *etap,
                                 Vector3f *wt) {
    float cosTheta_i = n.dot(wi);
    // Potentially flip interface orientation for Snell's law
    if (cosTheta_i < 0) {
        eta = 1 / eta;
        cosTheta_i = -cosTheta_i;
        n = -n;
    }

    // Compute $\cos\,\theta_\roman{t}$ using Snell's law
    float sin2Theta_i = std::max<float>(0.f, 1.f - (cosTheta_i*cosTheta_i));
    float sin2Theta_t = sin2Theta_i / (eta*eta);
    // Handle total internal reflection case
    if (sin2Theta_t >= 1)
        return false;

    float cosTheta_t = std::sqrt(1 - sin2Theta_t);

    *wt = -wi / eta + (cosTheta_i / eta - cosTheta_t) * Vector3f(n);
    // Provide relative IOR along ray to caller
    if (etap)
        *etap = eta;

    return true;
}



    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:

static Vector3f refract(const Vector3f& wi, const Vector3f& n, float eta)
{
    float cosThetaI = wi.dot(n);
    float sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
    float sin2ThetaT = eta * eta * sin2ThetaI;
    if (sin2ThetaT >= 1.0f) // full reflection
        return Vector3f(-wi[0], -wi[1], wi[2]);
    float cosThetaT = std::sqrt(1.0f - sin2ThetaT);
    if(cosThetaI > 0.f) cosThetaT = -cosThetaT;
    return eta * -wi + (eta * cosThetaI + cosThetaT) * n;
}
    float m_intIOR, m_extIOR;
    float m_eta, m_invEta;
    float m_transmittance;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END

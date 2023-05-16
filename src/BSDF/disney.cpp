/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DisneyBRDF : public BSDF {
public:
    DisneyBRDF(const PropertyList &propList)
    :baseColor (propList.getColor("baseColor")),
    roughness (propList.getFloat("roughness")),
    subsurface(propList.getFloat("subsurface")),
    metallic(propList.getFloat("metallic")),
    specular(propList.getFloat("specular")),
    specularTint(propList.getFloat("specularTint")),
    anystropic(propList.getFloat("anystropic")),
    sheen(propList.getFloat("sheen")),
    sheenTint(propList.getFloat("sheenTint")),
    clearcoat_(propList.getFloat("clearcoat")),
    clearcoatGloss(propList.getFloat("clearcoatGloss"))
    {}

    virtual bool isDelta()const{return false;}

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if(Frame::cosTheta(bRec.wi)<=0.f || Frame::cosTheta(bRec.wo)<=0.f)
        {
            return 0.f;
        }

        auto diffuse = DiffuseTerm(bRec.wi,bRec.wo,Normal3f(0.f,0.f,1.f));
        auto specular = SpecularTerm(bRec.wi,bRec.wo,Normal3f(0.f,0.f,1.f));
        auto clearcoat = ClearcoatTerm(bRec.wi,bRec.wo,Normal3f(0.f,0.f,1.f));
        return diffuse*(1.f-metallic) + specular + clearcoat;
    
    
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        //for disney brdf, all reflect, no transmit
        if(bRec.wo.z()<0.f || bRec.wi.z() < 0.f)
            return 0.f;
        
        float pdf_diffuse = PdfDiffuse(bRec.wi,bRec.wo);
        float pdf_specular = PdfSpecular(bRec.wi,bRec.wo);
        float pdf_clearcoat = PdfClearCoat(bRec.wi,bRec.wo);
        
        const float diffuse_prob = std::min(0.8f, 1.f - metallic);
        const float clearcoat_prob = clearcoat_ / (2.f+clearcoat_);
        
        return diffuse_prob*pdf_diffuse + (1.f-diffuse_prob)*( (1.f-clearcoat_prob)*pdf_specular + clearcoat_prob*pdf_clearcoat);
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord& bRec, Sampler* sampler) const
    {
        const float diffuse_prob = std::min(0.8f, 1.f - metallic);
        const float clearcoat_prob = clearcoat_ / (2.f+clearcoat_);
        if(sampler->next1D() < diffuse_prob)
        {//sample diffuse
            bRec.wo = SampleDiifuse(bRec.wi,sampler->next2D());
        }
        else{//sample specular or clearcoat
            if(sampler->next1D()> clearcoat_prob)
            {//specular
                bRec.wo = SampleSpecular(bRec.wi,sampler->next2D());
            }
            else{//clearcoat
                bRec.wo = SampleClearCoat(bRec.wi,sampler->next2D());
            }
        }

        float CosTheta = Frame::cosTheta(bRec.wo);
        if(CosTheta < 0.f)return 0.f;
        
        float pdf_ = pdf(bRec);

        if(pdf_ > 0.f)
            return eval(bRec)*CosTheta / pdf_;
    }






    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        
        float cosTheta = Frame::cosTheta(bRec.wi); 
        if (cosTheta <= 0.f)
            return Color3f(0.0f);
        bRec.measure = ESolidAngle;
        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */



        if(pdf(bRec) > 0.f)
            return eval(bRec)*cosTheta / pdf(bRec);
        else
            return 0.f;
    }






    bool isDiffuse() const {
        /* While disney BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }


    std::string toString() const {
        return tfm::format(
            "Disney[]"
        );
    }
private:
    float anystropic = 0.f; 
    float sheen = 0.f;
    float sheenTint = 0.f;
    float clearcoatGloss = 0.f;
    float clearcoat_ = 0.f;
    float roughness=0.f;
    float subsurface=0.f;
    float metallic=0.f;
    float specular=0.f;
    float specularTint=0.f;
    Color3f baseColor=0.f;
    //schlick's approximation for Fresnel
    static float SchlickFresnel(float CosTheta)
    {
        float tmp = std::clamp(1.f-CosTheta,0.f,1.f);
        float tmp2 = tmp*tmp;
        return tmp2*tmp2*tmp;
    }
    //helper function for DiffuseTerm
    std::pair<float,float> Diffuse_Parameter( const Vector3f& wi, const Vector3f& wo,  const Normal3f& n )const
    {
        Vector3f H = (wi+wo).normalized();
        const float wiDotH= H.dot(wi);
        const float Fd90 = 0.5f + 2.f * wiDotH * wiDotH * roughness;
        const float wiDotN = wi.z();//n.dot(wi);
        const float woDotN = wo.z();//n.dot(wo);
        const float Fi = SchlickFresnel(wiDotN);
        const float Fo = SchlickFresnel(woDotN);

        float Fd = (1.f + (Fd90-1.f)*Fi)*(1.f + (Fd90-1.f)*Fo);

        float Fss90 = wiDotH * wiDotH * roughness;
        float Fss = (1.f+(Fss90-1.f)*Fi)*(1.f+(Fss90-1.f)*Fo);
        float ss = 1.25f*(Fss*(1.f/(wiDotN+wiDotN)-0.5f)+0.5f);

        return {Fd,ss};
    }
    
    Color3f CorrectColor()const
    {
        Color3f Cdlin = baseColor;
        float Cdlum = 0.3f*Cdlin.r() + 0.6*Cdlin.g() +0.1*Cdlin.b();
        Color3f Ctint = (Cdlum>0.f)?(Cdlin/Cdlum):Color3f(1.f);
        Color3f Cspec = specular*( Color3f(1.f) + specularTint*(Ctint- Color3f(1.f)) ) ;
        Color3f Cspec0 = (0.08f*Cspec + metallic*(Cdlin-0.08f*Cspec) );
        return Cspec0;
    }

    static float GTR1(float NdotH, float alpha)
    {
        if(alpha>=1.f)return INV_PI;
        float alpha2=alpha*alpha;
        return (alpha2-1.f)/(M_PI*std::log(alpha2)*(1.f+(alpha2-1.f)*NdotH*NdotH) );
    }

    static float GTR2(float NdotH, float alpha)
    {
        float alpha2 = alpha*alpha;
        float tmp = 1.f+(alpha2-1.f)*NdotH*NdotH;
        return alpha2/(M_PI*tmp*tmp);
    }

    static float GTR2_ansio(float NdotH, float HdotX, float HdotY, float ax, float ay)
    {
        return 1.f/(M_PI*ax*ay*std::pow(std::pow(NdotH,2) + (std::pow(HdotX/ax,2)+std::pow(HdotY/ay,2)),2) );
    }

    static float Smith_GGX(float NdotWo,float alpha)
    {
        float a = alpha*alpha;
        float b = NdotWo*NdotWo;
        return 1.f/(NdotWo + std::sqrt(a+b-a*b));
    }

    static float Smith_GGX_aniso(float NdotWo,float Wo_x,float Wo_y, float alpha_x, float alpha_y)
    {
        return 1.f/(NdotWo + std::sqrt( std::pow(Wo_x*alpha_x, 2) + std::pow(Wo_y*alpha_y,2) + NdotWo*NdotWo ) );
    }
    
    //caculate DiffuseReflect term    
    Color3f DiffuseTerm( const Vector3f& wi, const Vector3f& wo,  const Normal3f& n )const
    {
        std::pair<float,float> res = Diffuse_Parameter(wi,wo,n);
        const float Fd = res.first;
        const  float ss = res.second;
        const Color3f Cdlin = baseColor;
        
        Color3f Cdlin = baseColor;
        float Cdlum = 0.3f*Cdlin.r() + 0.6*Cdlin.g() +0.1*Cdlin.b();
        Color3f Ctint = (Cdlum>0.f)?(Cdlin/Cdlum):Color3f(1.f);
        Color3f Csheen = (1.f-sheenTint)*Color3f(1.f)+sheenTint*Ctint;
        
        Vector3f H = (wi+wo).normalized();
        const float wiDotH = wi.dot(H);
        float FH = SchlickFresnel(wiDotH);
        Color3f Fsheen = FH*sheen*Csheen;

        Color3f diffuse =( Fd + (ss-Fd)*subsurface ) * Cdlin *INV_PI + Fsheen;
        return diffuse * (1.f - metallic);
    }
    //caculate SpecularReflect term
    Color3f SpecularTerm(const Vector3f& wi, const Vector3f& wo, const Normal3f& n)const
    {
        float alpha = roughness*roughness;
        float aspect = std::sqrt(1.f-0.9f*anystropic);
        float alpha_x = std::max(0.001f, alpha/aspect);
        float alpha_y = std::max(0.001f, alpha*aspect);
        
        const Vector3f H = (wi+wo).normalized();
        const float NdotH = H.z();//n.dot(H);
        const float woDotH = wo.dot(H);
        const float wiDotH = wi.dot(H);
        // float Ds = GTR2(NdotH,alpha);
        float Ds = GTR2_ansio(NdotH, H.x(), H.y(), alpha_x, alpha_y);

        float FH = SchlickFresnel(wiDotH);
        Color3f Cspec0 = CorrectColor();
        Color3f Fs = (1.f-FH) * Cspec0 + FH * Color3f(1.f);
        // float Gs = Smith_GGX(wi.z(),roughness)*Smith_GGX(wo.z(),roughness);
        // consider anisotropic effect
        float Gs = Smith_GGX_aniso(wi.z(),wi.x(),wi.y(),alpha_x,alpha_y)*Smith_GGX_aniso(wo.z(),wo.x(),wo.y(),alpha_x,alpha_y);

        Color3f specular = Gs * Fs * Ds ;// not divided by 4*wi.z()*wo.z() because the disney model don't need correcting

        return specular;
    }

    Color3f ClearcoatTerm(const Vector3f& wi, const Vector3f& wo, const Normal3f& n)const
    {
        const Vector3f H = (wi+wo).normalized();
        const float NdotH = H.z();//n.dot(H);
        float Dr = GTR1(NdotH, 0.1f*(1-clearcoatGloss) + 0.001f*clearcoatGloss );
        
        const float wiDotH = wi.dot(H);
        float FH = SchlickFresnel(wiDotH);

        float Fr = 0.04f*(1-FH)+1.f*FH;
        float WiDotN = wi.z();//n.dot(wi);
        float WoDotN = wo.z();//n.dot(wo);
        float Gr = Smith_GGX(WiDotN,0.25f) * Smith_GGX(WoDotN,0.25f);

        Color3f ClearcoatTerm = Gr * Fr * Dr * clearcoat_ * 0.25f;

        return ClearcoatTerm;

    }

    Vector3f SampleDiifuse(const Vector3f& wi, const Point2f& sample)const
    {
        float cosTheta = Frame::cosTheta(wi); 
        if (cosTheta <= 0.f)
            return Vector3f(0.f,0.f,0.f);
        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        Vector3f wo = Warp::squareToCosineHemisphere(sample);

        return wo;
    }

    float PdfDiffuse(const Vector3f& wi, const Vector3f& wo)const
    {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if ( Frame::cosTheta(wi) <= 0
            || Frame::cosTheta(wo) <= 0)
            return 0.0f;
        /* Importance sampling density wrt. solid angles:
           cos(theta) / pi.

           Note that the directions in 'bRec' are in local coordinates,
           so Frame::cosTheta() actually just returns the 'z' component.
        */
        return INV_PI * Frame::cosTheta(wo);
    }

    Vector3f SampleSpecular(const Vector3f& wi, const Point2f& sample)const
    {
        const float psi1 = sample[0];
        const float psi2 = sample[1];
        const float alpha = roughness*roughness;
        const float aspect = std::sqrt(1.f-0.9f*anystropic);
        const float alpha_x = std::max(0.001f, alpha/aspect);
        const float alpha_y = std::max(0.001f, alpha*aspect);
        
        float sin_phiH = alpha_y*std::sin(2.f*M_PI*psi1);
        float cos_phiH = alpha_x*std::cos(2.f*M_PI*psi1);
        const float R = std::sqrt(sin_phiH*sin_phiH+cos_phiH*cos_phiH);
        sin_phiH = sin_phiH/R;
        cos_phiH = cos_phiH/R;

        const float Omega_phiH = std::pow(sin_phiH/alpha_y,2) + std::pow(cos_phiH/alpha_x,2);
        const float cos_thetaH = std::sqrt( (Omega_phiH*(1.f-psi2)) /( (1.f-Omega_phiH)*psi2 + Omega_phiH));
        const float sin_thetaH = std::sqrt( std::max(0.f, 1.f-cos_thetaH*cos_thetaH));
        
        return Vector3f(cos_phiH*sin_thetaH, sin_phiH*sin_thetaH, cos_thetaH);
    
    }

    float PdfSpecular(const Vector3f& wi, const Vector3f& wo)const
    {
        const float alpha = roughness*roughness;
        const float aspect = std::sqrt(1.f-0.9f*anystropic);
        const float alpha_x = std::max(0.001f, alpha/aspect);
        const float alpha_y = std::max(0.001f, alpha*aspect);
        
        const Vector3f H = (wi+wo).normalized();
        const float NdotH = H.z();//n.dot(H);
        const float woDotH = wo.dot(H);
        const float wiDotH = wi.dot(H);
        // float Ds = GTR2(NdotH,alpha);
        float Ds = GTR2_ansio(NdotH, H.x(), H.y(), alpha_x, alpha_y);
        
        float pdf = Ds * NdotH / (4.f*wiDotH);
        return pdf;

    }

    Vector3f SampleClearCoat(const Vector3f& wi, const Point2f& sample) const
    {
        const float alpha = 0.1f*(1-clearcoatGloss) + 0.001f*clearcoatGloss;
        const float phi = 2.f*M_PI*sample[0];
        const float CosTheta = std::sqrt( (std::pow(alpha,2.f-2*sample[1]) / (alpha*alpha - 1.f) ) );
        const float SinTheta = std::sqrt( std::max(0.f,1.f-CosTheta*CosTheta) );

        return Vector3f(std::cos(phi)*SinTheta, std::sin(phi)*SinTheta, CosTheta);
    }

    float PdfClearCoat(const Vector3f& wi, const Vector3f& wo) const
    {
        Vector3f H = (wi+wo).normalized();
        float CosThetaH = H.z();
        float SinThetaH = std::sqrt( std::max(0.f,1.f-CosThetaH*CosThetaH) );
        float alpha = 0.1f*(1-clearcoatGloss) + 0.001f*clearcoatGloss;
        float alpha2 = alpha*alpha;
        float pdf = (alpha2-1.f)*CosThetaH*SinThetaH / ( (alpha2* CosThetaH*CosThetaH + SinThetaH*SinThetaH)* std::log(alpha) );

        return pdf;
    }   

};

NORI_REGISTER_CLASS(DisneyBRDF, "disney");
NORI_NAMESPACE_END

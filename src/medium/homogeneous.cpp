#include<nori/medium.h>
#include<nori/mesh.h>
#include"phasefunction.hpp"
#include<memory>

NORI_NAMESPACE_BEGIN

class Homogeneous:public Medium
{
    public:
    Homogeneous(const PropertyList& props)
    {
        float g = props.getFloat("g",0.0f);
        phase = std::make_shared<HenyneyGreestain>(g);
    }

    virtual Color3f transmitance(const Ray3f ray, Sampler *sampler) const override{}
    
    virtual Color3f sample(const Ray3f ray, Sampler *sampler, Intersection &mediumIts, Vector3f& newdir) const override{
        //sample channel
        int channel = 0;
        Vector3f pmf = SampleChannel(channel,sampler->next1D());
        //caculate sigma_t
        Color3f sigma_t = sigma_a + sigma_s;
        
        //sample distance
        float distance = -std::log(1.f - sampler->next1D())/sigma_t[channel];
        float t = std::min(distance, ray.maxt);
        
        mediumIts.p = ray(t);
        mediumIts.t = t;

        Color3f throughput=1.f;
        
        if(t > ray.maxt-1e-5f)
        { // transmit directly
            Color3f tr(std::exp(-t*sigma_t[0]),std::exp(-t*sigma_t[1]),std::exp(-t*sigma_t[2]));
            float p_surface = tr[0]*pmf[0]+tr[1]*pmf[1]+tr[2]*pmf[2];
            throughput = tr / p_surface;
            newdir = ray.d;
        }
        else{//scatter
            Color3f tr(std::exp(-t*sigma_t[0]),std::exp(-t*sigma_t[1]),std::exp(-t*sigma_t[2]));
            float p_scatter = tr[0]*pmf[0]*sigma_t[0]+tr[1]*pmf[1]*sigma_t[1]+tr[2]*pmf[2]*sigma_t[2];

            float pdf= phase->Sample_Direction(-ray.d,newdir,sampler->next2D());
            throughput = tr * sigma_s /p_scatter/pdf;
        }
    
        return throughput;


    }
    
    //interface to phase function
    virtual float SamplePhase(const Vector3f &wi, Vector3f& wo, const Point2f &sample) const override{
        return phase->Sample_Direction(wi, wo,sample);
    }
    //interface to phase function
    virtual float EvalPhase(const Vector3f &wi, const Vector3f &wo) const override{
        return phase->pdf(wi,wo);
    }
    
    virtual Color3f albedo() const override{}
    virtual Color3f get_sigma_s() const  override{return sigma_s;}
    
    EClassType getClassType() const { return EMedium; }

    private:
        std::shared_ptr<PhaseFunction> phase=nullptr;    
        Color3f sigma_a;
        Color3f sigma_s;

        Vector3f SampleChannel(int32_t& channel, const float sample)const
        {
            Color3f sigma_t = sigma_a+sigma_s;
            Color3f albedo = Color3f(sigma_s[0]/sigma_t[0],sigma_s[1]/sigma_t[1],sigma_s[2]/sigma_t[2]);

            Vector3f pmf = Vector3f(albedo[0],albedo[1],albedo[2]);
            pmf[0] = std::exp(pmf[0]);
            pmf[1] = std::exp(pmf[1]);
            pmf[2] = std::exp(pmf[2]);
            
            
            float sum = pmf[0]+pmf[1]+pmf[0];
            pmf /= sum;
            
            float psi = sample;
            if(psi>pmf[0]+pmf[1]){channel = 2;}
            else if(psi>pmf[0]){channel = 1;}
            else{channel = 0;}

            return pmf;
        }



};


NORI_NAMESPACE_END
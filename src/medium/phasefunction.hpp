#include<nori/common.h>
#include<nori/object.h>
#include<nori/vector.h>
#include<nori/sampler.h>
#include<nori/frame.h>

NORI_NAMESPACE_BEGIN
class PhaseFunction
{

    public:
    PhaseFunction() = default;
    ~PhaseFunction() = default;
    virtual float pdf(const Vector3f& vi, const Vector3f& vo)const = 0;
    virtual float Sample_Direction( const Vector3f& wi, Vector3f& wo,const Point2f& sample)const = 0;

};


class HenyneyGreestain: public PhaseFunction
{
    public:
    HenyneyGreestain() = delete;
    HenyneyGreestain(float g): g(g){}
    ~HenyneyGreestain()=default;
    

    virtual float pdf(const Vector3f& vi, const Vector3f& vo)const
    {
        //use HG phase function
        float cos_theta = vi.dot(vo)/(vi.norm()*vo.norm());
        float tmp = 1.f+g*g+2*g*cos_theta;
        float pdf = (1.f-g*g)/(4.f*M_PI*std::sqrt(tmp)*tmp);

        return pdf;

    }

    virtual float Sample_Direction(const Vector3f& wi, Vector3f& wo, const Point2f& sample)const
    {
        float cosTheta;
        if (std::abs(g) < 1e-3)
            cosTheta = 1 - 2 * sample[0];
        else {
            float sqrTerm = (1 - g * g) / (1 + g - 2 * g * sample[0]);
            cosTheta = -(1 + g * g - sqrTerm * sqrTerm) / (2 * g);
        }

        float sinTheta = std::sqrt(std::max((float)0, 1 - cosTheta * cosTheta));
        float phi = 2 * M_PI * sample[1];
        Frame f = Frame(wi);
        wo = f.toWorld(Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta));
        return pdf(wi,wo);
    }
    



    private:
        //asymetric parameter
        float g;
};

NORI_NAMESPACE_END
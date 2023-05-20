#pragma once

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/kdtree.h>

NORI_NAMESPACE_BEGIN

struct PhotonData
{
    uint8_t rgbe[4];
    uint8_t theta;
    uint8_t phi;

    static float m_cosTheta[0x100];
    static float m_sinTheta[0x100];
    static float m_cosPhi[0x100];
    static float m_sinPhi[0x100];
    static float m_expTable[0x100];
    static bool m_precomputeTableReady;



    PhotonData() = default;

    PhotonData(const Vector3f& dir, const Color3f& power)
    {
        if(!power.isValid())
        {
            std::cerr<<"Photon with invalid power:"<<power<<std::endl;
        }


        theta = uint8_t( std::clamp<int>( int(std::acos(dir.z()) * ( 256.f / M_PI ) ), 0, 255) );

        int tmp = std::min( int(std::atan2(dir.y(),dir.x())*( 256.f/(2.f*M_PI) )), 255 );
        if(tmp<0)
            phi = uint8_t(tmp+256);
        else
            phi = uint8_t(tmp);

        
        float max = power.maxCoeff();
        if(max<1e-32)
        {
            rgbe[0]=0;
            rgbe[1]=0;
            rgbe[2]=0;
            rgbe[3]=3;
        }
        else{
            int e=0;
            max=std::frexp(max,&e) * 256.f / max;

            rgbe[0] = uint8_t(power[0] * max);
            rgbe[1] = uint8_t(power[1] * max);
            rgbe[2] = uint8_t(power[2] * max);
            rgbe[3] = uint8_t(e+128);
        }

    }

    //tranform the power of photon to color
    inline Color3f getPower() const {
        return Color3f(rgbe[0],rgbe[1],rgbe[2]) * m_expTable[rgbe[3]];
    }

    //look up the direction of photon
    inline Vector3f getDirection() const
    {
        return Vector3f(m_cosPhi[phi] * m_sinTheta[theta],
                        m_sinPhi[phi] * m_sinTheta[theta],
                        m_cosTheta[theta]);

    }


    static bool intialize();

};





struct Photon: public GenericKDTreeNode<Point3f, PhotonData>
{
    public:
        Photon() {}
        explicit Photon(const Point3f& position, const Vector3f& dir, const Color3f& flux)
        :GenericKDTreeNode(position,PhotonData(dir,flux)) {}
        
        inline Vector3f getDirection()const{return data.getDirection();}
        
        inline Color3f getPower() const {return data.getPower();}

        inline Photon operator/ (float scale)const{
            return Photon(this->getPosition(), 
                        this->getDirection(), 
                        this->getPower() / scale);
        }


};

NORI_NAMESPACE_END
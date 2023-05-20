#include"photon.hpp"

NORI_NAMESPACE_BEGIN

bool PhotonData::intialize()
{
    for(int i=0;i<256;++i)
    {
        float angle = float(i) / 256.f * M_PI;
        m_cosTheta[i] = std::cos(angle);
        m_sinTheta[i] = std::sin(angle);
        m_cosPhi[i] = std::cos(2.f*angle);
        m_sinPhi[i] = std::sin(2.f*angle);
        m_expTable[i] = std::ldexp(1.f,i - (128+8) );
    }
    m_expTable[0]=0.f;
    m_precomputeTableReady=true;
    return true;
}



float PhotonData::m_cosTheta[0x100];
float PhotonData::m_sinTheta[0x100];
float PhotonData::m_cosPhi[0x100];
float PhotonData::m_sinPhi[0x100];
float PhotonData::m_expTable[0x100];

bool PhotonData::m_precomputeTableReady = PhotonData::intialize();






NORI_NAMESPACE_END
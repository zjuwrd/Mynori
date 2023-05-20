#include<nori/mesh.h>
#include<nori/transform.h>
#include<nori/common.h>
#include<nori/timer.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class Rectangle: public Mesh
{
public:
    Rectangle(const PropertyList& props){
        bool inverse_normal = props.getBoolean("inverse_normal",false);
        bool two_sided = props.getBoolean("two_sided",false);
        
        Timer timer;
        m_name = "Rectangle";
        Transform trafo = props.getTransform("toWorld");
     
        Point3f p0 = trafo*Point3f(-1.f,-1.f,0.f);
        Point3f p1 = trafo*Point3f(-1.f,1.f,0.f);
        Point3f p2 = trafo*Point3f(1.f,-1.f,0.f);
        Point3f p3 = trafo*Point3f(1.f,1.f,0.f);

        
        m_bbox.expandBy(p0);
        m_bbox.expandBy(p1);
        m_bbox.expandBy(p2);
        m_bbox.expandBy(p3);

        //fill m_V
        if(!two_sided)
            m_V.resize(3,4);
        else
            m_V.resize(3,8);
        m_V.col(0) = p3;
        m_V.col(1) = p2;
        m_V.col(2) = p0;
        m_V.col(3) = p1;

        

        //fill m_F
        if(!two_sided)
            m_F.resize(3,2);
        else
            m_F.resize(3,4);

        if(!inverse_normal)
        {
            m_F.col(0) = Vector3i( uint32_t(0),uint32_t(1),uint32_t(2)).cast<uint32_t>();
            m_F.col(1) = Vector3i(uint32_t(3),uint32_t(0),uint32_t(2)).cast<uint32_t>();
            
        }
        else{
            m_F.col(0) = Vector3i( uint32_t(2),uint32_t(1),uint32_t(0)).cast<uint32_t>();
            m_F.col(1) = Vector3i(uint32_t(2),uint32_t(0),uint32_t(3)).cast<uint32_t>();
        }


        if(two_sided)
        {
            Point3f q0 = trafo*Point3f(-1.f,-1.f,1e-5f);
            Point3f q1 = trafo*Point3f(-1.f,1.f,1e-5f);
            Point3f q2 = trafo*Point3f(1.f,-1.f,1e-5f);
            Point3f q3 = trafo*Point3f(1.f,1.f,1e-5f);

            m_bbox.expandBy(q0);
            m_bbox.expandBy(q1);
            m_bbox.expandBy(q2);
            m_bbox.expandBy(q3);
            
            m_V.col(4) = q3;
            m_V.col(5) = q2;
            m_V.col(6) = q0;
            m_V.col(7) = q1;

            
            if(!inverse_normal)
            {
                m_F.col(2) = Vector3i( uint32_t(6),uint32_t(5),uint32_t(4)).cast<uint32_t>();
                m_F.col(3) = Vector3i(uint32_t(6),uint32_t(4),uint32_t(7)).cast<uint32_t>();
                
            }
            else{
                m_F.col(2) = Vector3i( uint32_t(4),uint32_t(5),uint32_t(6)).cast<uint32_t>();
                m_F.col(3) = Vector3i(uint32_t(7),uint32_t(4),uint32_t(6)).cast<uint32_t>();
            }


        }





        cout <<m_name<<" loading..."
             << "done. (V=" << m_V.cols() << ", F=" << m_F.cols() << ", took "
             << timer.elapsedString() << " and "
             << memString(m_F.size() * sizeof(uint32_t) +
                          sizeof(float) * (m_V.size() + m_N.size() + m_UV.size()))
             << ")" << endl;
    

    }

};

NORI_REGISTER_CLASS(Rectangle, "rectangle")

NORI_NAMESPACE_END
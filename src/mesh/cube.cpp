#include<nori/mesh.h>
#include<nori/transform.h>
#include<nori/common.h>
#include<nori/timer.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class Cube: public Mesh
{
    public:
        Cube(const PropertyList& props )
        {
            bool inverse_normal = props.getBoolean("inverse_normal",false);
            Timer timer;
            m_name = "Cube";
            Transform trafo = props.getTransform("toWorld");
        
            const uint32_t m_face_count = 12;
            const uint32_t m_vertex_count = 24;

            m_V.resize(3,m_vertex_count);
            m_F.resize(3,m_face_count);
            m_UV.resize(2,m_vertex_count);


            std::vector<Point3f> vertices={
                {  1, -1, -1 }, {  1, -1,  1 }, { -1, -1,  1 }, { -1, -1, -1 },
                {  1,  1, -1 }, { -1,  1, -1 }, { -1,  1,  1 }, {  1,  1,  1 },
                {  1, -1, -1 }, {  1,  1, -1 }, {  1,  1,  1 }, {  1, -1,  1 },
                {  1, -1,  1 }, {  1,  1,  1 }, { -1,  1,  1 }, { -1, -1,  1 },
                { -1, -1,  1 }, { -1,  1,  1 }, { -1,  1, -1 }, { -1, -1, -1 },
                {  1,  1, -1 }, {  1, -1, -1 }, { -1, -1, -1 }, { -1,  1, -1 }    
            };


            std::vector<Normal3f> normals={
                { 0, -1,  0 }, {  0, -1,  0 }, {  0, -1,  0 }, {  0, -1,  0 }, {  0, 1, 0 },
                { 0,  1,  0 }, {  0,  1,  0 }, {  0,  1,  0 }, {  1,  0,  0 }, {  1, 0, 0 },
                { 1,  0,  0 }, {  1,  0,  0 }, {  0,  0,  1 }, {  0,  0,  1 }, {  0, 0, 1 },
                { 0,  0,  1 }, { -1,  0,  0 }, { -1,  0,  0 }, { -1,  0,  0 }, { -1, 0, 0 },
                { 0,  0, -1 }, {  0,  0, -1 }, {  0,  0, -1 }, {  0,  0, -1 }
            };
        
            std::vector<Point2f> texcoords={
                { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 },
                { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 },
                { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 },
                { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 }    
            };
        
            std::vector<Vector3i> triangles={
                {  0,  1,  2 }, {  3,  0,  2 }, {  4,  5,  6 }, {  7,  4,  6 },
                {  8,  9, 10 }, { 11,  8, 10 }, { 12, 13, 14 }, { 15, 12, 14 },
                { 16, 17, 18 }, { 19, 16, 18 }, { 20, 21, 22 }, { 23, 20, 22 }    
            };

            

            for(uint32_t i=0;i<m_vertex_count;++i)
            {
                //add vertex
                const auto& v_local = vertices[i];
                const auto& v_world = trafo * v_local;
                m_V.col(i) = v_world;

                //update the m_bbox
                m_bbox.expandBy(v_world);
                
                //add faces
                const auto& n_local = normals[i];
                const auto& n_world = (trafo * n_local).normalized();
                m_N.col(i) = inverse_normal ? -n_world : n_world;
                
                //add texcoords
                const auto& uv = texcoords[i];
                m_UV.col(i) = uv;

            }            


            cout <<m_name<<" loading..."
                << "done. (V=" << m_V.cols() << ", F=" << m_F.cols() << ", took "
                << timer.elapsedString() << " and "
                << memString(m_F.size() * sizeof(uint32_t) +
                            sizeof(float) * (m_V.size() + m_N.size() + m_UV.size()))
                << ")" << endl;
        


        }





};


NORI_NAMESPACE_END
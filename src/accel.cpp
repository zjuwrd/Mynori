/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
    
}

void Accel::build() {
    /* Nothing to do here for now */
    if(Use_Octree)
    {
        m_octree.clear();
        const uint32_t Primitive_count = m_mesh->getTriangleCount();
        size_t root_idx = CreateOctNode();
        for(uint32_t i=0;i<Primitive_count;++i)
        {
            Rootadd(i,root_idx);
        }
        adjust_octree();
    }
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    //Use Octree for acceleration
    if(Use_Octree)
    {
        foundIntersection = Octree_rayIntersect(ray_,its,f,shadowRay);
        if(shadowRay && foundIntersection)return true;
    }
    else
    {
        /* Brute force search through all triangles */
        for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
            float u, v, t;
            if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
                /* An intersection was found! Can terminate
                immediately if this is a shadow ray query */
                if (shadowRay)
                    return true;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                f = idx;
                foundIntersection = true;
            }
        }
    }
    
    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}


bool Accel::Octree_rayIntersect(const Ray3f &ray_, Intersection &its,uint32_t& primitive_idx_o, bool shadowRay, const size_t node_idx) const
{
    if(m_octree.empty())
    {
        return false;
    }
    const auto& node = m_octree[node_idx];


    Ray3f ray(ray_);
    bool foundIntersection = false;
    [[maybe_unused]]
    float farT=std::numeric_limits<float>::infinity();
    [[maybe_unused]]
    float nearT=-std::numeric_limits<float>::infinity();

    if(node.bbox.rayIntersect(ray))
    {
        if(node.children == -1)
        {
            for(const auto primitive_idx:node.primitives_idx)
            {
                float u, v, t;
                if (m_mesh->rayIntersect(primitive_idx,ray,u,v,t))
                {
                    if (shadowRay)
                        return true;
                    ray.maxt = its.t = t ;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    foundIntersection = true;
                    primitive_idx_o=primitive_idx;
                }
            }
        }
        else{
            for(int i=0;i<8;++i)
            {
                const auto& childnode = m_octree[node.children+i];
                if(childnode.bbox.rayIntersect(ray))
                {
                    if(Octree_rayIntersect(ray,its,primitive_idx_o,shadowRay,node.children+i))
                    {
                        if(shadowRay)return true;
                    }
                }
            }
        }
    }
    

}
    
    size_t Accel::CreateOctNode()
    {
        OctreeNode node;
        node.children = -1;
        m_octree.push_back(node);
        return m_octree.size() - 1;
    }

    size_t Accel::Create8OctNodes()
    {
        size_t idx = m_octree.size();
        for(int i = 0; i < 8; ++i)
        {
            CreateOctNode();
        }
        return idx;
    }

    size_t Accel::splitNode(size_t node_idx)
    {
        auto& node = m_octree[node_idx];
        node.children = Create8OctNodes();

        auto center = node.bbox.getCenter();
        auto min = node.bbox.min;
        auto max = node.bbox.max;
        
        for(int i=0;i<8;++i)
        {
            auto& children = m_octree[node.children+i];
            children.bbox = node.bbox;
            
            switch (i)
            {
            
            case 0:
                children.bbox.max = center;
                break;
            case 1:
                children.bbox.max.x() = center.x();
                children.bbox.max.y() = center.y();
                children.bbox.min.z() = center.z();
                break;
            case 2:
                children.bbox.max.x() = center.x();
                children.bbox.min.y() = center.y();
                children.bbox.max.z() = center.z();
                break;
            case 3:
                children.bbox.max.x() = center.x();
                children.bbox.min.y() = center.y();
                children.bbox.min.z() = center.z();
                break;
            case 4:
                children.bbox.min.x() = center.x();
                children.bbox.max.y() = center.y();
                children.bbox.max.z() = center.z();
                break;
            case 5:
                children.bbox.min.x() = center.x();
                children.bbox.max.y() = center.y();
                children.bbox.min.z() = center.z();
                break;
            case 6:
                children.bbox.min.x() = center.x();
                children.bbox.min.y() = center.y();
                children.bbox.max.z() = center.z();
                break;
            case 7:
                children.bbox.min.x() = center.x();
                children.bbox.min.y() = center.y();
                children.bbox.min.z() = center.z();
                break;
            default:
                break;
            }
            children.primitives_idx.reserve(node.primitives_idx.size());
            min = children.bbox.min;
            max = children.bbox.max;
            
            const auto& m_F = m_mesh->getIndices();
            const auto& m_V = m_mesh->getVertexPositions();

            for(auto primitive_idx : node.primitives_idx)
            {
                uint32_t i0 = m_F(0, primitive_idx), i1 = m_F(1, primitive_idx), i2 = m_F(2, primitive_idx);
                const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
                auto edge1 = p1 - p0;
                auto edge2 = p2 - p0;
                auto normal = edge1.cross(edge2);
                
                bool above= (Point3f(min[0],min[1],min[2]) - p0).dot(normal) > 0.f||(Point3f(max[0],min[1],min[2]) - p0).dot(normal) > 0.f||
                            (Point3f(min[0],max[1],min[2]) - p0).dot(normal) > 0.f||(Point3f(max[0],max[1],min[2]) - p0).dot(normal) > 0.f||
                            (Point3f(min[0],min[1],max[2]) - p0).dot(normal) > 0.f||(Point3f(max[0],min[1],max[2]) - p0).dot(normal) > 0.f||
                            (Point3f(min[0],max[1],max[2]) - p0).dot(normal) > 0.f||(Point3f(max[0],max[1],max[2]) - p0).dot(normal) > 0.f ;

                bool below = (Point3f(min[0],min[1],min[2]) - p0).dot(normal) <= 0.f||(Point3f(max[0],min[1],min[2]) - p0).dot(normal) <= 0.f||
                            (Point3f(min[0],max[1],min[2]) - p0).dot(normal) <= 0.f||(Point3f(max[0],max[1],min[2]) - p0).dot(normal) <= 0.f||
                            (Point3f(min[0],min[1],max[2]) - p0).dot(normal) <= 0.f||(Point3f(max[0],min[1],max[2]) - p0).dot(normal) <= 0.f||
                            (Point3f(min[0],max[1],max[2]) - p0).dot(normal) <= 0.f||(Point3f(max[0],max[1],max[2]) - p0).dot(normal) <= 0.f ;

                if(above&&below)
                {
                    children.primitives_idx.push_back(primitive_idx);
                }
            }
        }

        node.primitives_idx.clear();
        return node.children;
    }

    size_t Accel::Rootadd(size_t primitive_idx, size_t node_idx)
    {
        auto& node = m_octree[0];
        node.primitives_idx.push_back(primitive_idx);
        node.bbox.expandBy(m_mesh->getBoundingBox(primitive_idx));
    }

    void Accel::adjust_octree()
    {
        std::queue<size_t> adjustQue;
        if(m_octree.empty()||m_octree[0].primitives_idx.size()<=MaxPrimitivesPerNode)return;
        adjustQue.push(0);
        while(!adjustQue.empty())
        {
            const size_t node_idx = adjustQue.front();
            adjustQue.pop();
            auto& node = m_octree[node_idx];
            if(node.primitives_idx.size()<=MaxPrimitivesPerNode)
            {
                continue;
            }

            const size_t children_idx = splitNode(node_idx);
            for(size_t i=0;i<8;++i){adjustQue.push(children_idx+i);}
        }
    }



NORI_NAMESPACE_END


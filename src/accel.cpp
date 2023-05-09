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
        std::cout<<"building Octree ..."<<std::endl;
        m_octree.clear();
        m_octree.reserve(8*8*8*8*8);

        const uint32_t Primitive_count = m_mesh->getTriangleCount();
        size_t root_idx = CreateOctNode();

        for(uint32_t i=0;i<Primitive_count;++i)
        {
            Rootadd(i,root_idx);
        }
        std::cout<<"Root created."<<std::endl;
        adjust_octree();
        std::cout<<"Octree built."<<std::endl;

    }
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    //Use Octree for acceleration
    if(Use_Octree)
    {
        foundIntersection = Octree_rayIntersect(ray,its,f,shadowRay);
        if(shadowRay && foundIntersection)return true;
        // if(foundIntersection)std::cout<<"foundIntersect"<<std::endl;
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
        // std::cout<<"in"<<std::endl;
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


bool Accel::Octree_rayIntersect(Ray3f &ray_, Intersection &its,uint32_t& primitive_idx_o, bool shadowRay, const size_t node_idx) const
{
    if(m_octree.empty())
    {
        return false;
    }
    const auto& node = m_octree[node_idx];


    Ray3f& ray=ray_;
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
                        foundIntersection = true;
                        if(shadowRay)return true;
                    }
                }
            }
        }
    }
    
    return foundIntersection;

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
        // auto& node = m_octree[node_idx];
        size_t children_idx = Create8OctNodes();
        m_octree[node_idx].children = children_idx;

        auto center = m_octree[node_idx].bbox.getCenter();
        auto min = m_octree[node_idx].bbox.min;
        auto max = m_octree[node_idx].bbox.max;
        
        for(int i=0;i<8;++i)
        {
            auto& children = m_octree[m_octree[node_idx].children+i];
            children.bbox = m_octree[node_idx].bbox;
            const auto corner = children.bbox.getCorner(i);
            children.bbox.min[0] = std::min(corner[0], center[0]);
            children.bbox.min[1] = std::min(corner[1], center[1]);
            children.bbox.min[2] = std::min(corner[2], center[2]);
            children.bbox.max[0] = std::max(corner[0], center[0]);
            children.bbox.max[1] = std::max(corner[1], center[1]);
            children.bbox.max[2] = std::max(corner[2], center[2]); 


            children.primitives_idx.reserve(m_octree[node_idx].primitives_idx.size());
            min = children.bbox.min;
            max = children.bbox.max;
            
            [[maybe_unused]]
            const auto& m_F = m_mesh->getIndices();
            [[maybe_unused]]
            const auto& m_V = m_mesh->getVertexPositions();

            for(auto primitive_idx : m_octree[node_idx].primitives_idx)
            {
                auto prim_bbox = m_mesh->getBoundingBox(primitive_idx);
                if(prim_bbox.overlaps(children.bbox))
                {
                    children.primitives_idx.push_back(primitive_idx);
                }
                else{
                    continue;
                }


            }
        }

        m_octree[node_idx].primitives_idx.clear();
        return m_octree[node_idx].children;
    }

    size_t Accel::Rootadd(size_t primitive_idx, size_t node_idx)
    {
        auto& node = m_octree[0];
        node.primitives_idx.push_back(primitive_idx);
        node.bbox.expandBy(m_mesh->getBoundingBox(primitive_idx));
        return 0;

    }

    void Accel::adjust_octree()
    {
        std::queue<size_t> adjustQue;
        if(m_octree.empty()||m_octree[0].primitives_idx.size()<=MaxPrimitivesPerNode)return;
        adjustQue.push(0);
        size_t pivot = 0;
        int32_t depth = 0;
        while(!adjustQue.empty())
        {
            const size_t node_idx = adjustQue.front();
            adjustQue.pop();
            
    
            if(m_octree[node_idx].primitives_idx.size()<=MaxPrimitivesPerNode||depth>MaxOctreeDepth)
            {
                continue;
            }
            
            const size_t children_idx = splitNode(node_idx);
            for(size_t i=0;i<8;++i){
                
                adjustQue.push(children_idx+i);
            }
            
            if(pivot <= node_idx)
            {
                pivot = children_idx;
                depth++;
            }
        
        }
    }



NORI_NAMESPACE_END


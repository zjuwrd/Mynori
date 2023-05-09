/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/mesh.h>
#include<queue>

NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {

private:
    static constexpr int32_t MaxPrimitivesPerNode = 20;
    static constexpr int32_t MaxOctreeDepth = 8;

    const bool Use_Octree = true;

    struct OctreeNode {
        BoundingBox3f bbox;
        std::vector<uint32_t> primitives_idx;
        int32_t children;
    };
    
    std::vector<OctreeNode> m_octree;
    
    bool Octree_rayIntersect(Ray3f &ray_, Intersection &its,uint32_t& primitive_idx, bool shadowRay, const size_t node_idx = 0) const;
    
    size_t CreateOctNode();
    size_t Create8OctNodes();

    size_t splitNode(size_t node_idx);

    size_t Rootadd(size_t primitive_idx, size_t node_idx);

    void adjust_octree();


public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
};

NORI_NAMESPACE_END

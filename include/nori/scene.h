/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/accel.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Main scene data structure
 *
 * This class holds information on scene objects and is responsible for
 * coordinating rendering jobs. It also provides useful query routines that
 * are mostly used by the \ref Integrator implementations.
 */
class Scene : public NoriObject {
public:

    Color3f SampleLd(const Ray3f& ray, const Intersection& its, Sampler* sampler )const;


    ///return the number of emissive meshes
    size_t getEmissiveMeshesCount() const
    {
        size_t count = 0;
        for (auto mesh : m_meshes)
        {
            if (mesh->isEmitter())
                count++;
        }
        return count;
    }

    /// @brief return the idxth emissive mesh
    /// @param idx index of the emissive mesh
    /// @return idxth emissive mesh if it exists, else nullptr
    const Mesh* getLight(size_t idx)const
    {
        size_t count = 0;
        for (auto mesh : m_meshes)
        {
            if (mesh->isEmitter())
            {
                if (count == idx)
                    return mesh;
                count++;
            }
        }
        return nullptr;
    }

    std::vector<Mesh*> getLights()const
    {
        std::vector<Mesh*> ret;
        for(auto mesh : m_meshes)
        {
            if(mesh->isEmitter())ret.push_back(mesh);
        }

        return ret;
    }



    /// @brief Uniformly sample a light source
    /// @param point Uniform random number in [0,1]
    /// @param pdf possibility for the light source
    /// @return Light source sampled if exists, else nullptr
    const Mesh* SampleLight(const float point, float& pdf)const
    {
        size_t count = getEmissiveMeshesCount();
        if (count == 0)
            return nullptr;
        size_t idx = std::min((size_t)(point * count), count - 1);
        pdf = 1.0f / count;
        return getLight(idx);
    }

    const Mesh* SampleLight(const float point)const
    {
        size_t count = getEmissiveMeshesCount();
        if (count == 0)
            return nullptr;
        size_t idx = std::min((size_t)(point * count), count - 1);
        return getLight(idx);
    }




    /// Construct a new scene object
    Scene(const PropertyList &);

    /// Release all memory
    virtual ~Scene();

    /// Return a pointer to the scene's kd-tree
    const Accel *getAccel() const { return m_accel; }

    /// Return a pointer to the scene's integrator
    const Integrator *getIntegrator() const { return m_integrator; }

    /// Return a pointer to the scene's integrator
    Integrator *getIntegrator() { return m_integrator; }

    /// Return a pointer to the scene's camera
    const Camera *getCamera() const { return m_camera; }

    /// Return a pointer to the scene's sample generator (const version)
    const Sampler *getSampler() const { return m_sampler; }

    /// Return a pointer to the scene's sample generator
    Sampler *getSampler() { return m_sampler; }

    /// Return a reference to an array containing all meshes
    const std::vector<Mesh *> &getMeshes() const { return m_meshes; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its) const {
        return m_accel->rayIntersect(ray, its, false);
    }

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and \a only determine whether or not there is an intersection.
     *
     * This method much faster than the other ray tracing function,
     * but the performance comes at the cost of not providing any
     * additional information about the detected intersection
     * (not even its position).
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray) const {
        Intersection its; /* Unused */
        return m_accel->rayIntersect(ray, its, true);
    }

    /// \brief Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const {
        return m_accel->getBoundingBox();
    }

    /**
     * \brief Inherited from \ref NoriObject::activate()
     *
     * Initializes the internal data structures (kd-tree,
     * emitter sampling data structures, etc.)
     */
    void activate();

    /// Add a child object to the scene (meshes, integrators etc.)
    void addChild(NoriObject *obj);

    /// Return a string summary of the scene (for debugging purposes)
    std::string toString() const;

    EClassType getClassType() const { return EScene; }
private:
    std::vector<Mesh *> m_meshes;
    Integrator *m_integrator = nullptr;
    Sampler *m_sampler = nullptr;
    Camera *m_camera = nullptr;
    Accel *m_accel = nullptr;
};

NORI_NAMESPACE_END

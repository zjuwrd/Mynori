/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN




Mesh::Mesh() { }

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

std::pair<Point3f,Normal3f> Mesh::UniformSamplePoint(Sampler* sampler, float& pdf) const
{
    float s = sampler->next1D();
    const size_t f = m_dpdf->sample(s);
    pdf = 1.f/m_dpdf->getSum();


    const Point2f psi = sampler->next2D();
    const Point2f barycentric(1.f-std::sqrt(1.f-psi[0]),psi[1]*std::sqrt(1.f-psi[0]));

    uint32_t vi0 = m_F(0,f), vi1=m_F(1,f),vi2=m_F(2,f);
    const Point3f v0 = m_V.col(vi0), v1=m_V.col(vi1), v2=m_V.col(vi2);
    const Vector3f e1 = v1-v0;
    const Vector3f e2 = v2-v0;
    const Point3f pos = v0+e1*barycentric[0]+e2*barycentric[1];
    
    const Point3f normal = (m_V.cols() == m_N.cols())?(m_N.col(vi0)*(1.f-barycentric[0]-barycentric[1]) + m_N.col(vi1)*barycentric[1]+m_N.col(vi2)*barycentric[2]):
                                                        e1.cross(e2).normalized( );


    return {pos,normal};
}

void Mesh::samplePosition(const Point2f& sample, Point3f& position, Normal3f& normal, const float ratio)const
{
    float s = ratio;
    const size_t f = m_dpdf->sample(s);

    const Point2f psi = sample;
    const Point2f barycentric(1.f-std::sqrt(1.f-psi[0]),psi[1]*std::sqrt(1.f-psi[0]));

    uint32_t vi0 = m_F(0,f), vi1=m_F(1,f),vi2=m_F(2,f);
    const Point3f v0 = m_V.col(vi0), v1=m_V.col(vi1), v2=m_V.col(vi2);
    const Vector3f e1 = v1-v0;
    const Vector3f e2 = v2-v0;

    position = v0+e1*barycentric[0]+e2*barycentric[1];

    normal = (m_V.cols() == m_N.cols())?(m_N.col(vi0)*(1.f-barycentric[0]-barycentric[1]) + m_N.col(vi1)*barycentric[1]+m_N.col(vi2)*barycentric[2]):
                                                        e1.cross(e2).normalized( );


    
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }
    //Precomputing discrete pdf for sampling
    for(uint32_t f =0;f<m_F.cols();f++){
        float area = surfaceArea(f);
        m_dpdf->append(area);
    }


    m_dpdf->normalize();
    
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

MeshQeuryRecord Mesh::UniformSamplePoint(Sampler* sampler) const 
{
    MeshQeuryRecord result;
    uint32_t idx = m_dpdf->sample(sampler->next1D());
    Point2f rng = sampler->next2D();
    float alpha = 1 - sqrt(1 - rng.x());
    float beta = rng.y() * sqrt(1 - rng.x());
    Point3f v0 = m_V.col(m_F(0, idx));
    Point3f v1 = m_V.col(m_F(1, idx));
    Point3f v2 = m_V.col(m_F(2, idx));
    Point3f p = alpha * v0 + beta * v1 + (1 - alpha - beta) * v2;
    result.p = p;
    if (m_N.size() != 0) {
        Point3f n0 = m_N.col(m_F(0, idx));
        Point3f n1 = m_N.col(m_F(1, idx));
        Point3f n2 = m_N.col(m_F(2, idx));
        result.n = (alpha * n0 + beta * n1 + (1 - alpha - beta) * n2).normalized();
    } else {
        Vector3f e1 = v1 - v0;
        Vector3f e2 = v2 - v0;
        result.n = e1.cross(e2).normalized();
    }
    result.pdf = m_dpdf->getNormalization();
    return result;

}



NORI_NAMESPACE_END

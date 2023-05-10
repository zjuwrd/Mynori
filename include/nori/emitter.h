/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all emitters
 */

struct LightQueryRecord
{
    Point3f Light_Sample_point = Point3f(0.f,0.f,0.f);
    Point3f Primitive_point = Point3f(0.f,0.f,0.f);
    Normal3f AreaLight_normal = Point3f(0.f,1.f,0.f);
    
    LightQueryRecord(){}
    LightQueryRecord(Point3f Light_Sample_point,Point3f Primitive_point,Normal3f normal) 
    :Light_Sample_point(Light_Sample_point),Primitive_point(Primitive_point),AreaLight_normal(normal) {} 

};

class Emitter : public NoriObject {
public:
    virtual Color3f sample(const Mesh* mesh, LightQueryRecord& record, Sampler* sampler,float& pdf)const=0;
    virtual float pdf(const Mesh* mesh, LightQueryRecord& record)const=0;
    virtual Color3f getRadiance()const=0;
    virtual Color3f eval(const Mesh* mesh,const LightQueryRecord& record)const=0;
        


    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
};

NORI_NAMESPACE_END

#include<nori/emitter.h>
#include<nori/scene.h>

NORI_NAMESPACE_BEGIN



class AreaLight :public Emitter
{
    public:
        AreaLight(const PropertyList &props)
        :radiance(props.getColor("radiance",Color3f(1.0f))) {}
        
        /**
         * \brief Sample the emitter for a point on its surface
         * and return the incident radiance along the associated ray
         * */

        virtual Color3f sample(const Mesh* mesh, LightQueryRecord& record, Sampler* sampler, float& pdf)const
        {
            float light_pos_pdf=0.f;
            const auto result = mesh->UniformSamplePoint(sampler, light_pos_pdf);
            record.Light_Sample_point = result.first;
            record.AreaLight_normal = result.second;
            
            if (light_pos_pdf > 0.0f && !std::isnan(light_pos_pdf) && !std::isinf(light_pos_pdf)) {
                pdf = light_pos_pdf;
                return eval(mesh,record);            
            }

            pdf = 1.f;
            return 0.f;

        }
        virtual float pdf(const Mesh* mesh, LightQueryRecord& record)const
        {
        
            float cosTheta = record.AreaLight_normal.dot(record.Primitive_point-record.Light_Sample_point);
            //light to point
            if (cosTheta > 0.0f) {
                return mesh->getDpdf()->getNormalization();
            }

            return 0.0f;
        
        }

        virtual Color3f getRadiance()const{return radiance;}
        
        virtual Color3f eval(const Mesh* mesh, const LightQueryRecord& record)const
        {
            return record.AreaLight_normal.dot(record.Primitive_point - record.Light_Sample_point) >0.f?radiance:Color3f(0.f,0.f,0.f);

        }

        /**
         * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
         * provided by this instance
         * */
        std::string toString() const {
            return "AreaLight[]";
        }
    private:
        Color3f radiance;
};


NORI_REGISTER_CLASS(AreaLight, "area");

NORI_NAMESPACE_END
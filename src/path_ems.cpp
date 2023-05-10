#if 1
#include<nori/integrator.h>
#include<nori/scene.h>
#include<nori/emitter.h>
#include<nori/sampler.h>
#include<nori/warp.h>
#include<nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathEms : public Integrator {
    public:
    PathEms(const PropertyList &props) {}

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_)const{
    #if 1
        constexpr uint32_t MinDepth = 5;
        constexpr uint32_t MaxDepth = 100;
        Color3f throughput(1.0f);
        Color3f L(0.0f);
        Ray3f ray(ray_);
        uint32_t depth=1;
        Intersection its;
        while(depth<MaxDepth && scene->rayIntersect(ray,its))
        {
            if(its.mesh->isEmitter())
            {
                EmitterQueryRecord record(ray.o,its.p,its.shFrame.n);
                L+=throughput*its.mesh->getEmitter()->eval(record);
            }

            
            //sample directlight
            auto meshes = scene->getMeshes();
            Color3f Ld(0.f);

            for(const auto mesh:meshes)
            {
                if(mesh && mesh->isEmitter())
                {
                    if(mesh == its.mesh)
                        continue;
                    
                    EmitterQueryRecord record(its.p);
                    float pdf=0.f;
                    Color3f radiance=mesh->getEmitter()->sample(record,sampler);
                    pdf = record.pdf;
                    if(pdf>0.f)
                    {
                        Vector3f wi = record.wi;
                        Ray3f shadow_ray(its.p, wi, Epsilon, record.dist - Epsilon);
                        Intersection its_shadow;
                        if (!scene->rayIntersect(shadow_ray, its_shadow)) //no occlusion
                        {
                            BSDFQueryRecord bRec(its.shFrame.toLocal(wi),its.shFrame.toLocal(-ray.d),  EMeasure::ESolidAngle);
                            Ld += radiance * its.mesh->getBSDF()->eval(bRec);
                        }
                    }
                    
                }
            }
            L+=throughput*Ld;


            //Russian Roullete
            if(depth>MinDepth)
            {
                float p=std::min(throughput.maxCoeff(),0.99f);
                if(sampler->next1D() > p)
                    break;
                throughput/=p;
            }

            //sample indirectlight
            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
            Color3f bsdf=its.mesh->getBSDF()->sample(bRec,sampler->next2D());
            throughput*=bsdf;
            depth++;
            ray=Ray3f(its.p,its.shFrame.toWorld(bRec.wo));

        }
    
        return L/2.f;
    #endif
    
    

    }

    std::string toString() const {
        return "PathEms[]";
    }
};

NORI_REGISTER_CLASS(PathEms, "path_ems");

NORI_NAMESPACE_END
#endif


#if 0
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathEms : public Integrator {
 public:
  PathEms(const PropertyList& props) {}

  Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& _ray) const override {
    Color3f color = 0;//最终结果
    Color3f t = 1;//本次弹射对最终结果的贡献
    Ray3f rayRecursive = _ray;
    float probability;
    int depth = 1;
    int isDelta = 1;//判断是不是diffuse材质
    while (true) {
      Intersection its;
      if (!scene->rayIntersect(rayRecursive, its))
        break;
      if (its.mesh->isEmitter()) {//光源贡献
        LightQueryRecord lRecE(its.p,rayRecursive.o,  its.shFrame.n);
        color += t * its.mesh->getEmitter()->eval(its.mesh,lRecE) * isDelta;
      }
      if (its.mesh->getBSDF()->isDiffuse()) {
        float light_pdf;
        auto light = scene->SampleLight(sampler->next1D(),light_pdf);//均匀采样一个光源
        LightQueryRecord lRec;
        lRec.Light_Sample_point = its.p;

        float pdf;
        Color3f Li = light->getEmitter()->sample(light, lRec, sampler,pdf);//采样出射方向
        if (scene->rayIntersect(lRec.shadowRay())) {//出射方向有没有碰到物体
          Li = 0;
        }
        float cosTheta = Frame::cosTheta(its.shFrame.toLocal(lRec.wi()));
        BSDFQueryRecord bRec(its.toLocal(-rayRecursive.d), its.toLocal(lRec.wi()), ESolidAngle);
        Color3f f = its.mesh->getBSDF()->eval(bRec);//计算bsdf
        //color += Li * f * cosTheta * t;//叠加到最终结果上
        //2021.12.6更新：之前没有除以均匀随机选取光源的pdf
        color += Li * f * cosTheta * t / (1.0f / (float)scene->getEmissiveMeshesCount());
        isDelta = 0;
      } else {
        isDelta = 1;//不是diffuse就下次递归再说
      }
      //俄罗斯轮盘赌
      if (depth >= 3) {
        probability = std::min(t.maxCoeff(), 0.99f);
        if (sampler->next1D() > probability)
          break;
        t /= probability;
      }
      //采样下一个方向
      BSDFQueryRecord bRec(its.shFrame.toLocal(-rayRecursive.d));
      Color3f f = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
      t *= f;//叠加路径贡献
      rayRecursive = Ray3f(its.p, its.toWorld(bRec.wo));
      depth++;
    }
    return color;
  }

  std::string toString() const {
    return "PathEms[]";
  }
};

NORI_REGISTER_CLASS(PathEms, "path_ems");
NORI_NAMESPACE_END

#endif
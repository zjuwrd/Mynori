#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/sampler.h>
#include<nori/warp.h>
#include<nori/kdtree.h>
#include<nori/emitter.h>
#include<nori/timer.h>

#include"photon/photon.hpp"
#include<tbb/tbb.h>
#include<thread>
#include<memory>

#define USE_TBB

extern int threadCount;


NORI_NAMESPACE_BEGIN

class PM:public Integrator
{
    using PhotonMap = PointKDTree<Photon>;

    public:

    PM(const PropertyList &props)
    {
        m_photoncount = props.getInteger("photoncount",10000);
        m_radius = props.getFloat("radius",-1.f);
    }

    //shooting photons
    virtual void preprocess(const Scene* scene)
    {
        constexpr int Mindepth=3;
        constexpr int MaxDepth=100;

        Timer timer;    

        std::cout<<"Generating photons..."<<std::endl;
        std::unique_ptr<Sampler> Global_sampler(static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList())));

        m_photonmap = std::make_unique<PhotonMap>();
        m_photonmap->reserve(m_photoncount);


        if(m_radius<0.f)
        {
            m_radius = scene->getBoundingBox().getExtents().norm() / 500.f;
        }

        int shoot_cnt=0;

#ifdef USE_TBB
        tbb::mutex pm_mutex;
        
        std::thread PhotonThread(
        [&]{
                tbb::task_scheduler_init init(threadCount);
                shoot_cnt = tbb::parallel_reduce(tbb::blocked_range<size_t>(0,m_photoncount),0,
                [&](const tbb::blocked_range<size_t>& r,int init)->int
                {
                    
                    int local_shoot_count = 0;
                    const size_t local_photon_cnt = r.end()-r.begin();

                    std::unique_ptr<Sampler> sampler = Global_sampler->newClone();

                    // std::cout<<"thread "<<r.begin() <<" shooting "<<local_photon_cnt<<" photons"<<std::endl;


                    for(size_t photon_cnt= 0 ;photon_cnt < local_photon_cnt; )
                    {
                        const Mesh* emissiveMesh = scene->SampleLight(sampler->next1D());
                        const Emitter* emitter = emissiveMesh->getEmitter();

                        PhotonRay res = emitter->ShootPhoton(sampler.get());

                        if(res.success && !res.flux.isZero())
                        {
                            ++local_shoot_count;
                            Ray3f ray(res.ray);
                            Color3f flux(res.flux);
                            Intersection its;
                            int depth = 0;
                            Color3f throughput=1.f;
                        
                            while((depth<MaxDepth) && photon_cnt < local_photon_cnt && scene->rayIntersect(ray,its))
                            {
                                const BSDF* bsdf = its.mesh->getBSDF();
                                if(bsdf->isDiffuse())
                                {

                                    //mutex on m_phtonmap
                                    {
                                        tbb::mutex::scoped_lock lock(pm_mutex);
                                        m_photonmap->push_back(Photon(its.p, ray.d, throughput * flux));                                
                                        lock.release();
                                    }
                                    
                                    ++photon_cnt;
                                }

                                BSDFQueryRecord bQ(its.toLocal(-ray.d));

                                Color3f fr = bsdf->sample(bQ,sampler.get());
                                
                                Vector3f newdir = its.toWorld(bQ.wo);

                                throughput *= fr;

                                //Russian roullete
                                if(depth>Mindepth)
                                {
                                    float p = std::min(0.99f, throughput.maxCoeff());
                                    
                                    if(sampler->next1D() > p)
                                    {
                                        break;
                                    }                                
                                
                                    throughput /= p;
                                }   
                                ++depth;
                                ray = Ray3f(its.p,newdir,Epsilon,INFINITY);
                            }
                        }
                    }
                
                    return init+local_shoot_count;
                
                },[&](int x, int y)->int
                {
                    return x+y;
                });
            }
        );

    PhotonThread.join();
#else 
    const auto sampler = Global_sampler;
        for(int photon_cnt=0;photon_cnt<m_photoncount;)
        {
            const Mesh* emissiveMesh = scene->SampleLight(sampler->next1D());
            const Emitter* emitter = emissiveMesh->getEmitter();

            PhotonRay res = emitter->ShootPhoton(sampler);
            
            if(res.success && !res.flux.isZero())
            {
                ++shoot_cnt;
                Ray3f ray(res.ray);
                Color3f flux(res.flux);
                Intersection its;
                int depth = 0;

                Color3f throughput=1.f;


                while((depth<MaxDepth) && photon_cnt < m_photoncount && scene->rayIntersect(ray,its))
                {
                    const BSDF* bsdf = its.mesh->getBSDF();
                    if(bsdf->isDiffuse())
                    {
                        m_photonmap->push_back(Photon(its.p, ray.d, throughput * flux));
                        ++photon_cnt;
                    }

                    BSDFQueryRecord bQ(its.toLocal(-ray.d));
                    Color3f fr = bsdf->sample(bQ,sampler);

                    Vector3f newdir = its.toWorld(bQ.wo);

                    throughput *= fr;

                    //Russian roullete
                    if(depth>Mindepth)
                    {
                        float p = std::min(0.99f, throughput.maxCoeff());
                        if(sampler->next1D() > p)
                        {
                            break;
                        }
                        throughput /= p;
                    }
                    
                    ++depth;
                    ray = Ray3f(its.p,newdir,Epsilon,INFINITY);
                }

            }



        }

#endif

        std::cout<<"done. time elased "<<timer.elapsedString()<<"." <<std::endl;
        std::cout<<"collected photons / shooted photons = "<<m_photonmap->size()<<"/"<<shoot_cnt<<std::endl;

        timer.reset();
        std::cout<<"building photon map ..."<<std::endl;

        m_photonmap->scale(shoot_cnt);
        m_photonmap->build();

        std::cout<<"done. time elased "<<timer.elapsedString()<<"." <<std::endl;
        

    }


    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray_) const {
        constexpr int Mindepth = 3;
        constexpr int Maxdepth = 100;

        Ray3f curRay = ray_;
        int depth = 0;
        Intersection its;
        Color3f throughput = 1.f;
        Color3f L = 0.f;

        while(depth<Maxdepth && scene->rayIntersect(curRay,its) )
        {
            auto bsdf = its.mesh->getBSDF();
            
            if(its.mesh->isEmitter())
            {
                EmitterQueryRecord eQ(curRay.o, its.p, its.shFrame.n);
                const Emitter* emitter = its.mesh->getEmitter();
                L += throughput * emitter->eval(eQ);
            }


            if(bsdf->isDiffuse())
            {
                Color3f Lp = 0.f;

                std::vector<uint32_t> res;
                m_photonmap->search(its.p,m_radius, res);
                float area = M_PI * m_radius*m_radius;

                if(res.size() >0)
                {
                    for(const auto idx : res)
                    {
                        const auto& photon = (*m_photonmap)[idx];
                        BSDFQueryRecord bQ( its.toLocal(-curRay.d), its.toLocal(-photon.getDirection()) ,EMeasure::ESolidAngle);

                        Color3f fr = bsdf->eval(bQ);
                        float abscosTheta = std::abs( its.toLocal(photon.getDirection()).z());
                        
                        Lp += throughput *fr* abscosTheta * photon.getPower();
                        
                    }

                    Lp /= area;
                    L += Lp;
                }

                break;
            }

            //russian roullete
            if(depth>Mindepth)
            {
                float p = std::min(0.99f, throughput.maxCoeff());
                if(sampler->next1D()>p)
                {
                    break;
                }
                
                throughput /= p;
            }


            //sample new direction
            BSDFQueryRecord bQ(its.toLocal(-curRay.d));
            Color3f fr = bsdf->sample(bQ, sampler);
            Vector3f newdir = its.toWorld(bQ.wo);

            throughput *= fr;

            curRay = Ray3f(its.p, newdir);
            ++depth;

        }
    
        return L;
    
    }

    std::string toString() const {
        return "PMIntegrator[]";
    }

    private:
        float m_radius;
        int m_photoncount;

        std::unique_ptr<PhotonMap> m_photonmap=nullptr;


};




NORI_REGISTER_CLASS(PM, "PhotonMapping");


NORI_NAMESPACE_END

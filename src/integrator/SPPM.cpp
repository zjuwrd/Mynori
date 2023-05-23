#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/timer.h>
#include <nori/camera.h>
#include<nori/block.h>
#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <thread>
#include <memory>

#include "photon/photon.hpp"

#define USE_TBB

extern int threadCount;

NORI_NAMESPACE_BEGIN

class SPPM: public Integrator
{
    using PhotonMap = PointKDTree<Photon>;
    struct PixelMsg
    {
        Color3f flux;
        float radius;
        int photon_count;

        PixelMsg() = default;
        PixelMsg(float radius_):radius(radius_),photon_count(0), flux(0.f){}
        
    };



    public:
        SPPM(const PropertyList& props){
            m_photoncount = props.getInteger("PhotonCount",10000);
            m_iteration = props.getInteger("Iterations",10);
            alpha = props.getFloat("alpha",0.7);
            m_SharedRadius = props.getFloat("radius", 0.1f);
            m_emittedcount = 0;

        }

        virtual void preprocess(const Scene* scene)
        {
            m_photonmap = std::make_unique<PhotonMap>();
            m_photonmap->reserve(m_photoncount);            

            const Camera* camera = scene->getCamera();
            Vector2i outputSize = camera->getOutputSize();
            BlockGenerator blockGen(outputSize, NORI_BLOCK_SIZE);
            
            MsgImage.reserve(outputSize.y());

            for(int y=0;y<outputSize.y();++y)
            {
                MsgImage.push_back(std::vector<PixelMsg>(outputSize.x(),PixelMsg(m_SharedRadius)));
            }

        }
        
        virtual bool HasRenderMethod() const override {return true;}

        virtual void render(const Scene* scene,ImageBlock& image,  std::vector<NoriScreen*>& screens)override{

            for(int i=0;i<m_iteration; ++i)
            {
                image.clear();
                if(i%10 == 0)
                    std::cout<<"iteraion["<<i+1<<"]"<<" begin."<<std::endl; 
                
                ShootPass(scene);
                CapturePass(scene,image);
                if(i%10 == 0)
                    std::cout<<"iteration["<<i+1<<"]"<<" ends."<<std::endl;
            }

        }

        virtual Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray_) const override
        {
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
                    break;
                }


                if(bsdf->isDiffuse())
                {
                    Color3f Lp = 0.f;

                    std::vector<uint32_t> res;
                    m_photonmap->search(its.p,m_SharedRadius, res);
                    float area = M_PI * m_SharedRadius * m_SharedRadius;

                    if(res.size() >0)
                    {
                        for(const auto idx : res)
                        {
                            const auto& photon = (*m_photonmap)[idx];
                            BSDFQueryRecord bQ( its.toLocal(-curRay.d), its.toLocal(-photon.getDirection()) ,EMeasure::ESolidAngle);
                            Color3f fr = bsdf->eval(bQ);
                            float abscosTheta = std::abs( its.toLocal(photon.getDirection()).z());
                            Lp += throughput *fr* abscosTheta * photon.getPower() / m_emittedcount;   
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


        virtual Color3f estimate(const Scene* scene, Sampler* sampler, const Ray3f& ray_, const Vector2i pixel)  
        { 
            constexpr int Mindepth = 3;
            constexpr int Maxdepth = 100;

            Ray3f curRay = ray_;
            int depth = 0;
            Intersection its;
            Color3f throughput = 1.f;
            Color3f L = 0.f;

            auto& msg = MsgImage[pixel.y()][pixel.x()];

            while(depth<Maxdepth && scene->rayIntersect(curRay,its) )
            {
                auto bsdf = its.mesh->getBSDF();
                
                if(its.mesh->isEmitter())
                {
                    EmitterQueryRecord eQ(curRay.o, its.p, its.shFrame.n);
                    const Emitter* emitter = its.mesh->getEmitter();
                    L += throughput * emitter->eval(eQ);
                    break;
                }


                if(bsdf->isDiffuse())
                {
                    Color3f Lp = 0.f;

                    std::vector<uint32_t> res;
                    m_photonmap->search(its.p,m_SharedRadius, res);
                    float area = M_PI * m_SharedRadius * m_SharedRadius;

                    if(res.size() >0)
                    {
                        for(const auto idx : res)
                        {
                            const auto& photon = (*m_photonmap)[idx];
                            BSDFQueryRecord bQ( its.toLocal(-curRay.d), its.toLocal(-photon.getDirection()) ,EMeasure::ESolidAngle);
                            Color3f fr = bsdf->eval(bQ);
                            float abscosTheta = std::abs( its.toLocal(photon.getDirection()).z());
                            Lp += fr* abscosTheta * photon.getPower();   
                        }

                        float rate = (float)(msg.photon_count + alpha * res.size()) / (msg.photon_count + res.size());
                        msg.radius *= std::sqrt(rate);
                        msg.photon_count += res.size() *alpha;
                        msg.flux = (msg.flux + Lp) * rate ;
                    }

                    L += throughput * msg.flux / (M_PI * msg.radius * msg.radius * m_emittedcount);
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

        virtual std::string toString() const override
        {
            return tfm::format(
                "SPPM[\n"
                "]");
        }


    private:
        void CapturePass(const Scene* scene, ImageBlock& image){
            // defining Min-Max depth for Path tracing
            constexpr int Mindepth = 3;
            constexpr int Maxdepth = 100;
            // camra for the scene
            const Camera* camera = scene->getCamera();
            //get size of the image
            Vector2i outputSize = camera->getOutputSize();
            //Block Generator
            BlockGenerator blockGen(outputSize, NORI_BLOCK_SIZE);
            //Block range for tbb
            tbb::blocked_range<int> range(0, blockGen.getBlockCount());
            // parallel task for path tracing. Image is divided into blocks
            // to be dispatched to separate threads.

            auto map= [&](const tbb::blocked_range<int>& range)
            {
                //thread-local image block
                ImageBlock block(Vector2i(NORI_BLOCK_SIZE), camera->getReconstructionFilter());
                //thread-local sampler 
                std::unique_ptr<Sampler> sampler(scene->getSampler()->clone());

                for(int i=range.begin(); i < range.end(); ++i)
                {
                    //get next block to process
                    blockGen.next(block);
                    sampler->prepare(block);
                    
                    //process (progressively render) block
                    //get offset of the block
                    Point2i ofs = block.getOffset();
                    //get size of the block
                    Vector2i size = block.getSize();
                    block.clear();

                    //process the block by pixels
                    for(int y=0;y<size.y();++y)
                    {
                        for(int x=0;x<size.x();++x)
                        {
                            auto& msg = MsgImage[y][x];

                            //sample for each pixels couple of times predefined.
                            for(uint32_t samplecnt=0;samplecnt< 1/*sampler->getSampleCount()*/   ;++samplecnt)
                            {
                                //get positions in a single pixel
                                Point2f pixelSample = Point2f((float) (x + ofs.x()), (float) (y + ofs.y())) + sampler->next2D();
                                //get paerturesample
                                Point2f apertureSample = sampler->next2D();

                                /* Sample a ray from the camera */
                                Ray3f ray;
                                Color3f value = camera->sampleRay(ray, pixelSample, apertureSample);
                                /* Renew the Pixelmsgs and Compute the incident radiance */
#if 1
                                value *= estimate(scene,sampler.get(),ray,Vector2i(x+ ofs.x(),y+ofs.y()));
                                block.put(pixelSample, value);
#endif

#if 0
                                value *= Li(scene, sampler.get(), ray);
                                block.put(pixelSample, value);
#endif

                            }
                        
                        
                        
                        }
                    }

                    // result->put(block);
                    image.put(block);
                }
            };
            
            // do Capture in parallel
            tbb::task_scheduler_init init(threadCount);                
            tbb::parallel_for(range,map);

            // map(range);

        }
        
        void ShootPass(const Scene* scene){
            constexpr int Mindepth=3;
            constexpr int MaxDepth=100;

            Timer timer;    

            std::unique_ptr<Sampler> Global_sampler=scene->getSampler()->clone();

            m_photonmap->clear();
            m_photonmap->reserve(m_photoncount);

            int shoot_cnt=0;

#ifdef USE_TBB
            // mutex for photonmap
            tbb::mutex pm_mutex;
            //thread for shooting photons
            std::thread PhotonThread(
            [&]{
                    tbb::task_scheduler_init init(threadCount);
                    shoot_cnt = tbb::parallel_reduce(tbb::blocked_range<size_t>(0,m_photoncount),0,
                    [&](const tbb::blocked_range<size_t>& r,int init)->int
                    {
                        
                        int local_shoot_count = 0;
                        const size_t local_photon_cnt = r.end()-r.begin();

                        std::unique_ptr<Sampler> sampler = Global_sampler->newClone();
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

            timer.reset();
            // m_photonmap->scale(shoot_cnt);
            m_emittedcount += shoot_cnt;
            m_photonmap->build();

        }
        
        
        
        
        uint32_t m_photoncount;
        uint32_t m_emittedcount;
        uint32_t m_iteration;
        float m_SharedRadius;
        float alpha;
        std::unique_ptr<PhotonMap> m_photonmap=nullptr;
        // std::unique_ptr<ImageBlock> result = nullptr;

        std::vector< std::vector<PixelMsg> > MsgImage;
        

};

NORI_REGISTER_CLASS(SPPM,"sppm")

NORI_NAMESPACE_END
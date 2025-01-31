/*
    Stochastic Progressive photon mapping algorithm based on Nori
    implemented by wrd
*/

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
#include <nori/bitmap.h>
#include "photon/photon.hpp"

extern int threadCount;

NORI_NAMESPACE_BEGIN
struct PixelMsg
{
    Point2f pixel;
    float radius;
    Color3f flux;
    uint32_t p_nums;

    PixelMsg(){};

    PixelMsg(Point2f pixel, float radius, Color3f flux, uint32_t p_nums): 
        pixel(pixel), radius(radius), flux(flux), p_nums(p_nums) { }
};

class photon_sppm : public Integrator
{
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    
    photon_sppm(const PropertyList &props)
    {
        m_photonCount = props.getInteger("photonCount", (1<<14) );
        m_iteration = props.getInteger("iteration", 10 );
        m_sharedRadius = props.getFloat("radius", 0.2f);
        alpha = props.getFloat("alpha", 2.f/3.f);
        m_photonTotal = 0;

        record = props.getInteger("record",0);

    }

    virtual bool HasRenderMethod()const override{return true;}


    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray)const override{return 0.f;}

    virtual void render(const Scene *scene, ImageBlock &Image) override
    {
        
        // Initialize Pixel maps
        // size = Image.x * Image.y
        PixelMap.reserve(scene->getCamera()->getOutputSize().y());
        for(int y=0;y<scene->getCamera()->getOutputSize().y();++y)
        {
            PixelMap.push_back(std::vector< std::vector<PixelMsg> >( scene->getCamera()->getOutputSize().x() ) );
            for(int x=0;x<scene->getCamera()->getOutputSize().x();++x)
            {
                PixelMap[y].push_back(std::vector< PixelMsg>());
                for(int sp=0;sp<scene->getSampler()->getSampleCount();++sp)
                {
                    PixelMap[y][x].push_back(PixelMsg(Point2f(x,y),m_sharedRadius,0.f,0));
                }
            }    
        }

        // print out messages        
        std::cout
            << "\nsample nums: " << PixelMap.size()*PixelMap[0].size()*PixelMap[0][0].size()
            << "\niteration nums: " << m_iteration
            << "\nphoton nums per pass: " << m_photonCount
            << std::endl;
        
        
        // generate sampler
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
        
        const Camera *camera = scene->getCamera();

        for (uint32_t i = 0; i < m_iteration; ++i)
        {
            Timer timer;
            m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
            m_photonMap->reserve(m_photonCount); 

            uint32_t storedPhotons = 0;
            uint32_t photonEmitter = 0;
/*-----------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------Photon pass--------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------------------------*/
            std::unique_ptr<Sampler> Global_sampler=scene->getSampler()->clone();

            m_photonMap->clear();
            m_photonMap->reserve(m_photonCount);

            int shoot_cnt=0;

            // mutex for photonmap
            tbb::mutex pm_mutex;
            //thread for shooting photons
            std::thread PhotonThread(
            [&]{
                    // accelerate photon pass via tbb
                    tbb::task_scheduler_init init(threadCount);
                    shoot_cnt = tbb::parallel_reduce(tbb::blocked_range<size_t>(0,m_photonCount),0,
                    [&](const tbb::blocked_range<size_t>& r,int init)->int
                    {
                        
                        int local_shoot_count = 0;
                        const size_t local_photon_cnt = r.end()-r.begin();

                        std::unique_ptr<Sampler> sampler = Global_sampler->newClone();
                        for(size_t photon_cnt= 0 ;photon_cnt < local_photon_cnt; )
                        {
                            // Sample light
                            const Mesh* emissiveMesh = scene->SampleLight(sampler->next1D());
                            const Emitter* emitter = emissiveMesh->getEmitter();
                            // Shooting photons with emitter
                            PhotonRay res = emitter->ShootPhoton(sampler.get());
                            // photons being emitted successfully, update the photon map
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
                                            m_photonMap->push_back(Photon(its.p, ray.d, throughput * flux));                                
                                            lock.release();
                                        }
                                        ++photon_cnt;
                                    }

                                    // sample new directions via BSDF
                                    BSDFQueryRecord bQ(its.toLocal(-ray.d));
                                    Color3f fr = bsdf->sample(bQ,sampler.get());
                                    Vector3f newdir = its.toWorld(bQ.wo);
                                    // update throughput
                                    throughput *= fr;
                                    //Russian roullete
                                    if(depth>MinDepth)
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

        timer.reset();
        // m_photonmap->scale(shoot_cnt);
        m_photonTotal += shoot_cnt;
        m_photonMap->build();

        std::cout<<"Photon Map pass done."<<std::endl;

/*-----------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------PT pass------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------------------------*/
// defining Min-Max depth for Path tracing
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
                            for(int sp=0;sp<scene->getSampler()->getSampleCount();++sp)
                            {
                                //sample for each pixels couple of times predefined.
                                {
                                    //get positions in a single pixel
                                    Point2f pixelSample = Point2f((float) (x + ofs.x()), (float) (y + ofs.y())) + sampler->next2D();
                                    //get paerturesample
                                    Point2f apertureSample = sampler->next2D();

                                    /* Sample a ray from the camera */
                                    Ray3f ray;
                                    Color3f value = camera->sampleRay(ray, pixelSample, apertureSample);
                                    /* Renew the Pixelmsgs and Compute the incident radiance */
                                    
                                    {
                                        //get the pixelmsg for the pixel
                                        auto& pmsg = PixelMap[y+ofs.y()][x+ofs.x()][sp];
                                        //sample point inside on pixel
                                        Point2f sample = pmsg.pixel + sampler->next2D();

                                        Intersection its;
                                        Color3f throughput(1.f);
                                        uint32_t depth = 0;

                                        
                                        for ( ;depth<MaxDepth && scene->rayIntersect(ray,its);++depth)
                                        {
                                            //Le
                                            if (its.mesh->isEmitter())
                                            {
                                                EmitterQueryRecord eRec(ray.o, its.p, its.shFrame.n);
                                                Color3f power = its.mesh->getEmitter()->eval(eRec);
                                                block.put(pixelSample, power * throughput);
                                                break;
                                            }

                                            //diffusive surface
                                            if (its.mesh->getBSDF()->isDiffuse())
                                            {
                                                // Query the photon map
                                                std::vector<uint32_t> QueryRes;
                                                m_photonMap->search(its.p, pmsg.radius, QueryRes);
                                                //No photons
                                                if (QueryRes.size() == 0)
                                                    break;
                                                //Compute the new radius and flux
                                                float rate = (float)(pmsg.p_nums + alpha * QueryRes.size()) / (pmsg.p_nums + QueryRes.size());
                                                pmsg.p_nums += QueryRes.size() * alpha;
                                                pmsg.radius = pmsg.radius * std::sqrt( rate);
                                                Color3f Lp(0.f);
                                                // Compute the radiance
                                                for (auto idx : QueryRes)
                                                {
                                                    Photon &photon = (*m_photonMap)[idx];
                                                    BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(-photon.getDirection()), ESolidAngle);
                                                    Color3f fr = its.mesh->getBSDF()->eval(bRec);
                                                    float abscosTheta = std::max(0.f, bRec.wo.z()) ;
                                                    fr *= abscosTheta;
                                                    Lp += fr * photon.getPower();
                                                }
                                                // Update the flux message
                                                pmsg.flux = (pmsg.flux + Lp * throughput) * rate;
                                                break;
                                            }

                                            //specular surface
                                            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
                                            // sample next direction using BSDF
                                            Color3f fr = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                                            if (fr.maxCoeff() == 0.f)
                                                break;
                                            // update throughput
                                            throughput *= fr;
                                            Ray3f ro(its.p, its.shFrame.toWorld(bRec.wo));
                                            memcpy((void*)&ray,(void*)&ro,sizeof(Ray3f));
                                            
                                            // Russian roulette
                                            if (depth > MinDepth)
                                            {
                                                float q = throughput.maxCoeff();
                                                if (sampler->next1D() > q)
                                                    break;
                                                throughput /= q;
                                            }
                                        
                                        }
                                        Color3f power = pmsg.flux / (m_photonTotal * M_PI * pmsg.radius * pmsg.radius);
                                        value *= power;
                                    }
                                                                        
                                    block.put(pixelSample, value);

                                }
                        
                            }
                        
                        }
                    }

                    // result->put(block);
                    Image.put(block);
                }
            };
            
            // do Capture in parallel
            tbb::task_scheduler_init init(threadCount);                
            tbb::parallel_for(range,map);


            // if(i%10==9)
            // {
            //     cout << "(the "<<i-8<<"-"<< i+1 <<" passes took " << timer.elapsedString() << ")" << endl;
            // }
            if(record)
            {
                std::unique_ptr<Bitmap> bitmap(Image.toBitmap());
                std::string outputName = "results/sppm/SSPM" + std::to_string(i+1) ;
                bitmap->savePNG(outputName);
            }

        }//iter结束
    }






    virtual std::string toString() const override
    {
        return tfm::format(
            "StochasticProgressivePhotonMapper[\n"
            "]");
    }

private:
    static constexpr int MinDepth=3;
    static constexpr int MaxDepth=100;

    bool record = false;
    uint32_t m_photonCount;                   
    uint32_t m_photonTotal;                   
    uint32_t m_iteration;                     
    float m_sharedRadius;                     
    float alpha;                            
    std::vector< std::vector< std::vector<PixelMsg> > >  PixelMap; 
    std::unique_ptr<PhotonMap> m_photonMap; 
};

NORI_REGISTER_CLASS(photon_sppm, "sppm");
NORI_NAMESPACE_END
/*
    Debiased Progressive photon mapping algorithm based on Nori
    implemented by wrd
*/
#include <nori/integrator.h>
#include <nori/scene.h>
#include<nori/sampler.h>
#include<nori/warp.h>
#include<nori/kdtree.h>
#include<nori/emitter.h>
#include<nori/timer.h>
#include<nori/camera.h>
#include<nori/block.h>
#include <nori/gui.h>

#include"photon/photon.hpp"
#include<tbb/tbb.h>
#include<thread>
#include<memory>

#define DEBIASED

extern int threadCount;

NORI_NAMESPACE_BEGIN

class DebiasedPPM: public Integrator
{
    struct EstimateRes
    {
        Color3f Li0 = 0.f;
        Color3f Li1= 0.f;
        Color3f Lik= 0.f;
    };


    using PhotonMap = PointKDTree<Photon>;

    public:
        DebiasedPPM(const PropertyList &props)
        {
            init_k = props.getInteger("initk",1);
            init_Rad = props.getFloat("radius",0.2f);
            iterations = props.getInteger("iterations",1);
            PhotonsUnit = props.getInteger("photonsunit",100000);

        }

        virtual bool HasRenderMethod()const override{return true;}


        virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray_) const override {
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


        virtual void render(const Scene* scene, ImageBlock& image, std::vector<NoriScreen*>& screens) override
        {
            std::unique_ptr<Sampler> sampler = scene->getSampler()->newClone();
            
            KImg = std::make_unique<ImageBlock>(scene->getCamera()->getOutputSize(),scene->getCamera()->getReconstructionFilter());
            KImg->clear();
            NoriScreen* Kscreen = new NoriScreen(*KImg, "Biased Result");
            screens.push_back(Kscreen);

            for(int i=0;i<iterations;++i)
            {
                float alpha_pmf =0.f;
                const uint32_t k = init_k+i;
                const uint32_t j0 = SampleJ(k, sampler->next1D(), alpha_pmf);
                const uint32_t j1 = j0+1;

                std::cout<<"------------------------------------------------"<<std::endl;
                std::cout<<"iteration "<<i+1<<" starts."<<std::endl;
                std::cout<<"j="<<j0<<std::endl;
                std::cout<<"k="<<k<<std::endl;


                rad0 = init_Rad * std::pow(float(j0),(alpha-1.f)/2.f);
                rad1 = init_Rad * std::pow(float(j1),(alpha-1.f)/2.f);
                radk = init_Rad * std::pow(float(k),(alpha-1.f)/2.f);
                std::cout<<"\nrad0="<<rad0<<std::endl;
                std::cout<<"rad1="<<rad1<<std::endl;
                std::cout<<"radk="<<radk<<std::endl;
                
                std::cout<<"Photon Pass"<<std::endl;
                PhotonPass(j0, scene);
                std::cout<<"PT Pass"<<std::endl;
                CapturePass(i,j0,alpha_pmf,scene,image, screens);
                std::cout<<"done."<<std::endl;
                std::cout<<"totoal photon count: "<<total_photons<<std::endl;
            }
        }

        virtual std::string toString()const{
            return "DebiasedPPM[]";
        }


    private:
        void CapturePass(const uint32_t iteration, const uint32_t j, const float pmf, const Scene* scene, ImageBlock& image,std::vector<NoriScreen*>& screens){
            // defining Min-Max depth for Path tracing
            constexpr int Mindepth = 3;
            constexpr int Maxdepth = 50;
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

            BlockGenerator KblockGen(outputSize, NORI_BLOCK_SIZE);
            
            m_radius = init_Rad;

            auto map= [&](const tbb::blocked_range<int>& range)
            {
                //thread-local image block
                ImageBlock block(Vector2i(NORI_BLOCK_SIZE), camera->getReconstructionFilter());
                ImageBlock Kblock(Vector2i(NORI_BLOCK_SIZE), camera->getReconstructionFilter());    

                //thread-local sampler 
                std::unique_ptr<Sampler> sampler(scene->getSampler()->clone());

                for(int i=range.begin(); i < range.end(); ++i)
                {
                    //get next block to process
                    blockGen.next(block);
                    KblockGen.next(Kblock);
                    // PosDeltablockGen.next(PosDeltablock);
                    // NegDeltablockGen.next(NegDeltablock);

                    sampler->prepare(block);
                    
                    //process (progressively render) block
                    //get offset of the block
                    Point2i ofs = block.getOffset();
                    //get size of the block
                    Vector2i size = block.getSize();
                    block.clear();
                    Kblock.clear();
                    
                    //process the block by pixels
                    for(int y=0;y<size.y();++y)
                    {
                        for(int x=0;x<size.x();++x)
                        {
                            //sample for each pixels couple of times predefined.
                            for(uint32_t samplecnt=0;samplecnt< sampler->getSampleCount();++samplecnt)
                            {
                                //get positions in a single pixel
                                Point2f pixelSample = Point2f((float) (x + ofs.x()), (float) (y + ofs.y())) + sampler->next2D();
                                //get paerturesample
                                Point2f apertureSample = sampler->next2D();

                                /* Sample a ray from the camera */
                                Ray3f ray;
                                Color3f value = camera->sampleRay(ray, pixelSample, apertureSample);
                                /* Renew the Pixelmsgs and Compute the incident radiance */
                                
                            #ifdef DEBIASED
                                auto res = estimate(j,pmf, scene, sampler.get(),ray);
                                //final image
                                Color3f ImageVal = value * (res.Lik + (res.Li1 - res.Li0) / pmf)/5.3f;
                                ImageVal = ImageVal.clamp();
                                //biased k image
                                Color3f KimageVal = value * res.Lik;
                                KimageVal = KimageVal.clamp();
                                //delta graph
                                block.put(pixelSample, ImageVal );
                                Kblock.put(pixelSample, res.Lik);

                                // if(DeltaimageVal.minCoeff() > 0)
                                //     PosDeltablock.put(pixelSample, DeltaimageVal);
                                // else if(DeltaimageVal.maxCoeff() < 0)
                                //     NegDeltablock.put(pixelSample, -DeltaimageVal);

                            
                            #else
                                value *= Li(scene,sampler.get(),ray);
                                block.put(pixelSample ,value);
                            #endif
                            }
                        }
                    }
                #ifdef DEBIASED
                    // image.weighted_put(block, 1.f/( float(iteration) + 1.f) );
                    image.put(block);
                    KImg->put(Kblock);
                    // PosDeltaImg->put(PosDeltablock);
                    // NegDeltaImg->put(NegDeltablock);
                #else
                    image.put(block);
                #endif
                }
            };
            
            // do Capture in parallel
            tbb::task_scheduler_init init(threadCount);                
            tbb::parallel_for(range,map);

            // map(range);

        }



        EstimateRes estimate(const uint32_t j, const float pmf, const Scene *scene, Sampler *sampler, const Ray3f &ray_) const {
            
            EstimateRes res;

            Ray3f curRay = ray_;
            int depth = 0;
            Intersection its;
            Color3f throughput = 1.f;

            while(depth<MaxDepth && scene->rayIntersect(curRay,its) )
            {
                auto bsdf = its.mesh->getBSDF();
                
                if(its.mesh->isEmitter())
                {
                    EmitterQueryRecord eQ(curRay.o, its.p, its.shFrame.n);
                    const Emitter* emitter = its.mesh->getEmitter();
                    res.Lik += throughput * emitter->eval(eQ);
                    break;
                }


                if(bsdf->isDiffuse())
                {
                    Color3f Lp0 = 0.f, Lp1=0.f, Lpk=0.f;

                    std::vector<uint32_t> res0;
                    std::vector<uint32_t> res1;
                    std::vector<uint32_t> resk;
                    
                    //estimate j0
                    m_photonmap->search(its.p,rad0, res0);
                    float area0 = M_PI * rad0 * rad0;
                    if(res0.size() >0)
                    {
                        for(const auto idx : res0)
                        {
                            const auto& photon = (*m_photonmap)[idx];
                            BSDFQueryRecord bQ( its.toLocal(-curRay.d), its.toLocal(-photon.getDirection()) ,EMeasure::ESolidAngle);

                            Color3f fr = bsdf->eval(bQ);
                            float abscosTheta = std::abs( its.toLocal(photon.getDirection()).z());
                            Lp0 += throughput *fr* abscosTheta * photon.getPower();

                        }
                        Lp0 /= area0;
                        res.Li0 += Lp0;
                    }

                    //estimate j1
                    m_photonmap->search(its.p,rad1, res1);
                    float area1 = M_PI * rad1 * rad1;
                    if(res1.size() >0)
                    {
                        for(const auto idx : res1)
                        {
                            const auto& photon = (*m_photonmap)[idx];
                            BSDFQueryRecord bQ( its.toLocal(-curRay.d), its.toLocal(-photon.getDirection()) ,EMeasure::ESolidAngle);
                            Color3f fr = bsdf->eval(bQ);
                            float abscosTheta = std::abs( its.toLocal(photon.getDirection()).z());
                            Lp1 += throughput *fr* abscosTheta * photon.getPower();
                        }
                        Lp1 /= area1;
                        res.Li1 += Lp1;
                    }

                    //estimate k
                    m_photonmap->search(its.p,radk, resk);
                    float areak = M_PI * radk * radk;
                    if(resk.size() >0)
                    {
                        for(const auto idx : resk)
                        {
                            const auto& photon = (*m_photonmap)[idx];
                            BSDFQueryRecord bQ( its.toLocal(-curRay.d), its.toLocal(-photon.getDirection()) ,EMeasure::ESolidAngle);
                            Color3f fr = bsdf->eval(bQ);
                            float abscosTheta = std::abs( its.toLocal(photon.getDirection()).z());
                            Lpk += throughput *fr* abscosTheta * photon.getPower();
                        }
                        Lpk /= areak;
                        res.Lik += Lpk;
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
        
            // return res.Lik;
            return res;

            

        }



        void PhotonPass(const uint32_t j, const Scene* scene)
        {
            uint32_t PhotonCount = (std::pow( (j+1), 1-c*alpha))  * PhotonsUnit;
            total_photons += PhotonCount;
            //control the size of photon map
            if(PhotonCount > 80000000) PhotonCount = 80000000;
            Timer timer;
            std::unique_ptr<Sampler> Global_sampler(static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList())));
            m_photonmap = std::make_unique<PhotonMap>();
            m_photonmap->reserve(PhotonCount);

            int shoot_cnt=0;

            tbb::mutex pm_mutex;            
            std::thread PhotonThread(
            [&]{
                    tbb::task_scheduler_init init(threadCount);
                    shoot_cnt = tbb::parallel_reduce(tbb::blocked_range<size_t>(0,PhotonCount),0,
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
        timer.reset();
        m_photonmap->scale(shoot_cnt);
        m_photonmap->build();
        }

        uint32_t SampleJ(const uint32_t k, const float psi, float& pmf) const
        {
            const float k_regular = std::pow(k+0.f,alpha-1.f);
            float n = k*std::pow((1.f-psi),1.f/(alpha-1.f));

            uint32_t n1 = uint32_t(std::max(float(k) , std::floor(n)));
            uint32_t n2 = n1+1;

            float alpha_pmf = (std::pow( float(n1),alpha-1.f) - std::pow( float(n2),alpha-1.f))/k_regular;
            if(alpha_pmf<=0.f)
            {
                alpha_pmf =1.f;
            }

            pmf = alpha_pmf;

            return n1;
        
        }


        uint32_t SampleJ(const float psi, float& pmf )const
        {
            const float k_regular = std::pow(init_k+0.f,alpha-1.f);
            float n = init_k*std::pow((1.f-psi),1.f/(alpha-1.f));

            uint32_t n1 = uint32_t(std::max(float(init_k) , std::floor(n)));
            uint32_t n2 = n1+1;

            float alpha_pmf = (std::pow( float(n1),alpha-1.f) - std::pow( float(n2),alpha-1.f))/k_regular;
            if(alpha_pmf<=0.f)
            {
                alpha_pmf =1.f;
            }

            pmf = alpha_pmf;

            return n1;
        }



        //init k for Prime Estimizer
        float rad0,rad1,radk;

        uint32_t total_photons=0;
        uint32_t init_k=0;
        uint32_t iterations=0;
        float init_Rad =  0.2f;
        float m_radius = 0.2f;
        uint32_t PhotonsUnit = 100000;
        std::unique_ptr<PhotonMap> m_photonmap=nullptr;
        std::unique_ptr<ImageBlock> KImg=nullptr;
        std::unique_ptr<ImageBlock> PosDeltaImg=nullptr;
        std::unique_ptr<ImageBlock> NegDeltaImg=nullptr;
        static constexpr int Mindepth=3;
        static constexpr int MaxDepth=100;
        static constexpr float alpha = 2.f/3.f;
        static constexpr float c = 1.00001f;


};



NORI_REGISTER_CLASS(DebiasedPPM, "debiasedppm")

NORI_NAMESPACE_END
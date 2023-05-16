#pragma once

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

#include<queue>
#include "photon.hpp"


NORI_NAMESPACE_BEGIN
class PhotonMap
{
    private:
    class KDTree
    {
        private:
            struct Node{
                public:
                int32_t axis;
                int32_t photon_idx;  //idx for photon
                int32_t left_ch=INVALID_IDX;//default to -1 meaning null
                int32_t right_ch=INVALID_IDX;//default to -1 meaning null
                int32_t depth=0;
                Node(int32_t axis){this->axis = axis;}
                Node(int32_t axis, int32_t photon_idx, int32_t left_ch, int32_t right_ch)
                :photon_idx(photon_idx),left_ch(left_ch),right_ch(right_ch)
                {}
                ~Node()=default;
            };

            static constexpr int32_t INVALID_IDX = -1;
            std::vector<Node> nodes={}; // first element is empty
            std::vector<int32_t> photon_indices;
            const Photon* photons;
            const int32_t photon_count;
            

            inline int32_t create_node(int32_t depth)
            {
                nodes.push_back(Node(depth%3));
                return nodes.size()-1;    
            }

            int32_t build_KDTree(int32_t* part_photon_indices,int32_t nphotons, int32_t depth);

        private:

            void KNeareast(const Point3f& point, int32_t search_idx, std::vector<int32_t>& maxheap, const size_t K) const;

            void RangeSearch(const Point3f& point, int32_t search_idx, std::queue<int32_t>& list, float radius) const;

        public:
            KDTree(const Photon* photons, const int32_t photon_count)
            :photons(photons),photon_count(photon_count){}
            
            std::vector<int32_t> find_K_Nearest(const Point3f& point, const size_t K)const;
            std::queue<int32_t> Sphere_Search(const Point3f& point, const float radius)const;
            void construct_KDTree();
    };

    public:
        PhotonMap(){}

        inline void SetPhotons(const std::vector<Photon>& photons){ this->photons = photons; }
        inline void add(const Photon& photon){ photons.push_back(photon); }
        inline bool MapEmpty(){ return photons.empty();}

        void BuildMap();
        void ReleaseMap();

        const std::vector<Photon>& get_Photons() const { return photons; }
        
        std::vector<int32_t> find_K_Nearest_Photons(const Point3f& point, const size_t K) const;
        std::vector<int32_t> find_Photons_in_Sphere(const Point3f& point, const float radius) const;
    
    private:
        std::vector<Photon> photons;
        std::unique_ptr<KDTree> kdtree=nullptr;

};
    

NORI_NAMESPACE_END
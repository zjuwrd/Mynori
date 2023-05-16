#include"photonmap.hpp"

#include<vector>
#include<stack>
#include<algorithm>

NORI_NAMESPACE_BEGIN

int32_t PhotonMap::KDTree::build_KDTree(int32_t* part_photon_indices,int32_t nphotons, int32_t depth)
{
    if(nphotons <= 0)return INVALID_IDX;
    int axis = depth%3;
    auto criteria = [&](const int32_t id1, const int32_t id2){return photons[id1].position[axis] < photons[id2].position[axis]; };
    std::sort(part_photon_indices,part_photon_indices+nphotons, criteria );

    int32_t mid_idx = part_photon_indices[nphotons/2];
    int32_t new_node_idx = create_node(depth);
    auto& new_node = nodes[new_node_idx];
    new_node.photon_idx = mid_idx;
    
    new_node.left_ch = build_KDTree(part_photon_indices,nphotons/2,depth+1);
    new_node.right_ch = build_KDTree(part_photon_indices+nphotons/2,nphotons-nphotons/2-1,depth+1);

    return new_node_idx;
}

void PhotonMap::KDTree::KNeareast(const Point3f& point, int32_t search_idx, std::vector<int32_t>& maxheap, const size_t K)const
{
    auto MaxHeapCriteria = ([&](int32_t idx1, int32_t idx2){ 
        const auto& photon1 = photons[idx1];
        const auto& photon2 = photons[idx2];
        float distance1_quadr = (photon1.position - point).dot(photon1.position - point);
        float distance2_quadr = (photon2.position - point).dot(photon2.position - point);    
        return distance1_quadr < distance2_quadr;
    });
    
    if(search_idx<0)return;

    const Node& node = nodes[search_idx];
    const Photon& photon = photons[node.photon_idx];

    maxheap.push_back(node.photon_idx);
    std::push_heap(maxheap.begin(),maxheap.end(),MaxHeapCriteria);

    if(maxheap.size()>K)
    {
        std::pop_heap(maxheap.begin(),maxheap.end(),MaxHeapCriteria);
        maxheap.pop_back();
    }

    auto axis = node.axis;
    bool choose_left = point[axis] < photon.position[axis];

    if(choose_left)
    {
        KNeareast(point, node.left_ch, maxheap,K);
    }            
    else{
        KNeareast(point, node.left_ch, maxheap,K);    
    }

    const Photon& Kth_photon = photons[maxheap[0]];
    if( (point-Kth_photon.position).norm() > std::abs(point[axis]-photon.position[axis]) )
    {
        if(choose_left)
        {
            KNeareast(point, node.right_ch,maxheap,K);
        }
        else{
            KNeareast(point, node.left_ch,maxheap,K);
        }
    }
}

void PhotonMap::KDTree::RangeSearch(const Point3f& point, int32_t search_idx, std::queue<int32_t>& list, float radius)const
{
    if(search_idx<0)return;

    const Node& node = nodes[search_idx];
    const Photon& photon = photons[node.photon_idx];
    
    const Point3f delta = photon.position-point;
    if(delta.dot(delta) > radius*radius)
    {
        list.push(node.photon_idx);
    }
    
    float dist = point[node.axis]- photon.position[node.axis];
    bool choose_left = dist<0.f;
    if(choose_left)
    {
        RangeSearch(point,node.left_ch,list,radius);
    }
    else{
        RangeSearch(point,node.right_ch,list,radius);
    }

    if(dist*dist>radius*radius)
    {
        if(choose_left)
        {
            RangeSearch(point,node.right_ch,list,radius);
        }
        else{
            RangeSearch(point,node.left_ch,list,radius);
        }
    
    }

}

std::vector<int32_t> PhotonMap::KDTree::find_K_Nearest(const Point3f& point, const size_t K) const
{
    std::vector<int32_t> maxheap;
    KNeareast(point,0,maxheap,K);
    return maxheap;
}

std::queue<int32_t> PhotonMap::KDTree::Sphere_Search(const Point3f& point, const float radius)const
{
    std::queue<int32_t> list;
    RangeSearch(point,0,list,radius);
    return list;
}

void PhotonMap::KDTree::construct_KDTree()
{
    photon_indices.resize(photon_count);
    for(int32_t i=0;i<photon_count;++i)
    {
        photon_indices[i]=i;
    }
    build_KDTree(photon_indices.data(),photon_count,0);
}


void  PhotonMap::BuildMap(){
    kdtree = std::make_unique<KDTree>(KDTree(photons.data(),photons.size()));
    kdtree->construct_KDTree();   
}

void  PhotonMap::ReleaseMap(){
    kdtree.reset();
    photons.clear();
}

std::vector<int32_t>  PhotonMap::find_K_Nearest_Photons(const Point3f& point, const size_t K)const
{
    return kdtree->find_K_Nearest(point ,K);
}

std::vector<int32_t>  PhotonMap::find_Photons_in_Sphere(const Point3f& point, const float radius)const
{
    std::queue que=kdtree->Sphere_Search(point,radius);
    std::vector<int32_t> ret;
    while(!que.empty())
    {
        ret.push_back(que.front());
        que.pop();
    }

    return ret;
}

NORI_NAMESPACE_END
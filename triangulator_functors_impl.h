/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef TRIANGULATOR_FUNCTORS_IMPL_H
#define	TRIANGULATOR_FUNCTORS_IMPL_H

#include <algorithm>
#include <unordered_set>
#include "triangulator.h"
#include "triangulator_types_impl.h"
#include "utilities.h"


template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::simplex_t_hash {
    public:
        std::size_t operator()(const simplex_t& s) const {
            return vector_hash<edge_t,edge_t_hash<skeleton_t>>()(s.edges);
        }
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_t_comparison {
    public:
        bool operator()(const edge_t e, const edge_t f) const {
            auto min_e = std::min(e.m_target,e.m_source);
            auto min_f = std::min(f.m_target,f.m_source);
            if(min_e < min_f)
                return true;
            if(min_e > min_f)
                return false;
            return std::max(e.m_target,e.m_source) < std::max(f.m_target,f.m_source);
        }
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_t_equality {
    private:
        const edge_t e;
    public:
        edge_t_equality(const edge_t edge) : e(edge) { }
        bool operator()(const edge_t f) const {
            auto min_e = std::min(e.m_target,e.m_source);
            auto min_f = std::min(f.m_target,f.m_source);
            auto max_e = std::max(e.m_target,e.m_source);
            auto max_f = std::max(f.m_target,f.m_source);
            return (min_e == min_f) && (max_e == max_f);
        }
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_subsets_generator {
    private:
        typedef edge_container_t triangle_faces_t;
        typedef std::pair<edge_container_t,std::vector<triangle_faces_t>> tetrahedron_faces_t;
        typedef std::pair<edge_container_t,std::vector<tetrahedron_faces_t>> pentachoron_faces_t;             
        std::vector<tetrahedron_faces_t> tetrahedra;
        
        edge_container_t generateSubsets(const vertex_container_t& vertices, const skeleton_t& skeleton) {
            if(vertices.size() == 2) {
                auto e_pair = boost::edge(*vertices.begin(),*boost::next(vertices.begin()),skeleton);
                if(!e_pair.second)
                    throw std::logic_error("Generated non-existent edge in simplex subsets. Probably due to a faulty hypersurface triangulation.");
                edge_container_t edge;
                edge.push_back(e_pair.first);
                return edge;
            }
            edge_container_t edges;
            std::vector<edge_container_t> temp_triangles;
            for(auto iter = vertices.begin(); iter != vertices.end(); ++iter) {
                vertex_container_t temp_vertices;
                temp_vertices.push_back(*iter);
                auto next = boost::next(iter);
                for(int i = 1; i < vertices.size()-1; ++i) {
                    if(next == vertices.end())
                        next = vertices.begin();
                    temp_vertices.push_back(*next);
                    ++next;
                }
                std::sort(temp_vertices.begin(),temp_vertices.end());
                edge_container_t temp_edges = this->generateSubsets(temp_vertices,skeleton);
                if(temp_edges.size() == 3) {
                    temp_triangles.push_back(temp_edges);
                }
                for(auto in_iter = temp_edges.begin(); in_iter != temp_edges.end(); ++in_iter) {
                    if(std::find_if(edges.begin(), edges.end(), edge_t_equality(*in_iter)) == edges.end()) {
                        edges.push_back(*in_iter);
                    }                    
                }
            }        
            if(edges.size() == 6) {
                tetrahedra.push_back(std::make_pair(edges,temp_triangles));
            }
            return edges;
        }
        void printFaces(const pentachoron_faces_t& penta) {
            typedef edge_container_t triangle_faces_t;
            typedef std::pair<edge_container_t,std::vector<triangle_faces_t>> tetrahedron_faces_t;
            typedef std::pair<edge_container_t,std::vector<tetrahedron_faces_t>> pentachoron_faces_t;  
            std::cout << "Pentachoron: (";
            for(auto iter = penta.first.begin(); iter != penta.first.end(); ++iter) {
                std::cout << *iter;
                if(std::next(iter) != penta.first.end())
                    std::cout << ",";
            }
            std::cout << ")" << std::endl;
            for(auto iter = penta.second.begin(); iter != penta.second.end(); ++iter) {
                std::cout << " Tetrahedron: (";
                for(auto in_iter = iter->first.begin(); in_iter != iter->first.end(); ++in_iter) {
                    std::cout << *in_iter;
                    if(std::next(in_iter) != iter->first.end())
                        std::cout << ",";
                }
                std::cout << ")" << std::endl;
                for(auto in_iter = iter->second.begin(); in_iter != iter->second.end(); ++in_iter) {
                    std::cout << "  Triangle: (";
                    for(auto tri_iter = in_iter->begin(); tri_iter != in_iter->end(); ++tri_iter) {
                        std::cout << *tri_iter;
                        if(std::next(tri_iter) != in_iter->end())
                            std::cout << ",";
                    }
                    std::cout << ")" << std::endl;
                }
            }
            
        }
    public:
        edge_subsets_generator() {
            
        }        
        pentachoron_faces_t operator()(const vertex_container_t& vertices, const skeleton_t& skeleton) {
            edge_container_t penta_edges = this->generateSubsets(vertices,skeleton);
            pentachoron_faces_t penta_pair(penta_edges,tetrahedra);
            return penta_pair;
            
        }
};

#endif	/* TRIANGULATOR_FUNCTORS_IMPL_H */


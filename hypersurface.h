/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef HYPERSURFACE_H
#define	HYPERSURFACE_H

#include <boost/graph/adjacency_list.hpp>
#include "triangulator.h"
#include "utilities.h"

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::hypersurface_t::hypersurface_t(tetrahedron_container_t& vec) {
    for(auto iter = vec.begin(); iter != vec.end(); ++iter) {
                this->addTetrahedron(*iter);
    }
    vec.clear();
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::hypersurface_t::clear() {
    triangle_set.clear();
    tetrahedron_set.clear();
    hypersurface_skeleton.clear();
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::hypersurface_t::addTetrahedron(const vertex_t v1, const vertex_t v2, const vertex_t v3, const vertex_t v4) {
    vertex_container_t vec;
    vec.push_back(v1);
    vec.push_back(v2);
    vec.push_back(v3);
    vec.push_back(v4);
    this->addTetrahedron(vec);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::hypersurface_t::addTetrahedron(vertex_container_t& vec) {
    if(vec.size() != 4)
        throw std::logic_error("Trying to add a tetrahedron with != 4 vertices.");

    //sort and check internal vertex index
    std::sort(vec.begin(),vec.end());
    int i = 0;
    for(auto iter = vec.begin(); iter != vec.end(); ++iter) {
        auto index_iter = vertex_index_bimap.left.find(*iter);
        if(index_iter != vertex_index_bimap.left.end())
            *iter = index_iter->second;
        else {
            vertex_t internal_index = boost::add_vertex(hypersurface_skeleton);
            vertex_index_bimap.insert(vertex_indices_t(*iter,internal_index));
            *(iter) = internal_index;
        }
    }

    //check if tetrahedron is already added
    if(tetrahedron_set.find(vec) != tetrahedron_set.end())
        return ; 

    //add tetrahedron to set, triangles to set and edges to graph
    tetrahedron_set.insert(vec);
    for(auto iter = vec.begin(); iter != vec.end(); ++iter) {
        for(auto in_iter = boost::next(iter); in_iter != vec.end(); ++in_iter) {
            if(!boost::edge(*iter,*in_iter,hypersurface_skeleton).second)
                boost::add_edge(*iter,*in_iter,hypersurface_skeleton);
        }

        vertex_container_t triangle_vec;
        triangle_vec.push_back(*iter);
        for(auto in_iter = boost::next(iter); triangle_vec.size() < 3; ++in_iter) {
            if(in_iter==vec.end())
                in_iter = vec.begin();
            triangle_vec.push_back(*in_iter);
        }
        std::sort(triangle_vec.begin(),triangle_vec.end());
        triangle_set.insert(triangle_vec);
    }
} 

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::hypersurface_t::printSimplices() {
    std::cout << "Triangles: " << std::endl;
    for(auto iter = triangle_set.begin(); iter != triangle_set.end(); ++iter) {
        std::cout << "(";
        for(auto in_iter = iter->begin(); in_iter != iter->end(); ++in_iter) {
            std::cout << vertex_index_bimap.right.find(*in_iter)->second;
            if(std::next(in_iter) != iter->end())
                std::cout << ",";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << "Tetrahedra: " << std::endl;
    for(auto iter = tetrahedron_set.begin(); iter != tetrahedron_set.end(); ++iter) {
        std::cout << "(";
        for(auto in_iter = iter->begin(); in_iter != iter->end(); ++in_iter) {
            std::cout << vertex_index_bimap.right.find(*in_iter)->second;
            if(std::next(in_iter) != iter->end())
                std::cout << ",";
        }
        std::cout << ")" << std::endl;
    }
}

#endif	/* HYPERSURFACE_H */


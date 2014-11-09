/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef TRIANGULATOR_IMPL_H
#define	TRIANGULATOR_IMPL_H

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/graph/smallest_last_ordering.hpp>
#include "hypersurface.h"
#include "triangulator.h"
#include "utilities.h"


template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
template<typename map_t>
typename map_t::mapped_type& triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getOrCreateSimplex(map_t& map, const edge_container_t& edges) {
    typename map_t::iterator iter = map.emplace(edges,edges).first;
    return (*iter).second;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
template<typename type_t> 
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::addSimplexToContainer(std::vector<std::reference_wrapper<type_t> >& container, type_t& type) {
    if(std::find(container.begin(),container.end(),type) == container.end())
        container.push_back(std::ref(type));
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::triangulator_t(hypersurface_t& hypersurface, const geodesics_t& geodesics) : vertex_index_bimap(hypersurface.vertex_index_bimap) {
    null_vertex = skeleton_descriptors::null_vertex();
    if(!this->checkHypersurface(hypersurface))
        throw std::logic_error("Check of hypersurface failed. Maybe it has a boundary?");
    boost::copy_graph(hypersurface.hypersurface_skeleton,this->skeleton);
    this->copySimplices(hypersurface);
    hypersurface.clear();
    
    this->generateTriangulation();
    this->setInitialLengths(geodesics);
    
    vertex_hyper_index_pred<skeleton_t> h_pred(2,skeleton);
    auto filtered_skeleton = this->getFilteredSkeleton(h_pred);
    for(auto iter = boost::vertices(filtered_skeleton).first; iter != boost::vertices(filtered_skeleton).second; ++iter) {
        skeleton[*iter].evolved = false;
    }
}
  
template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_t 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getVertexAt(hyper_index_t h, vertex_t v) const {
    vertex_t temp_v = v;
    hyper_index_t actual_h = this->getVertexProperty(v).hyper_index;
    
    if(h > actual_h) {
        for(int i = 0; i < (h-actual_h); ++i) 
            temp_v = this->getVertexProperty(temp_v).next;
    }
    else {
        for(int i = 0; i < (actual_h-h); ++i) 
            temp_v = this->getVertexProperty(temp_v).last;
    }
    return temp_v;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setInitialLengths(const geodesics_t& geodesics) {
    for(auto iter = boost::edges(skeleton).first; iter != boost::edges(skeleton).second; ++iter) {
        vertex_t v0 = this->getSource(*iter);
        vertex_t v1 = this->getTarget(*iter);
        
        
        const hyper_index_t h0 = this->getVertexProperty(v0).hyper_index;
        const hyper_index_t h1 = this->getVertexProperty(v1).hyper_index;
        if(std::abs(h0 - h1) > 1)
            throw std::logic_error("Trying to set edge length between Sigma0 and Sigma2.");
        
        this->setEdgeLengthSq(*iter,geodesics(h0,vertex_index_bimap.right.find(this->getVertexAt(0,v0))->second,h1,vertex_index_bimap.right.find(this->getVertexAt(0,v1))->second));
    }
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::~triangulator_t() {
    //destroy everything
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
bool triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::checkHypersurface(const hypersurface_t& hypersurface) const {
    typedef typename hypersurface_t::traits::vertex_iterator v_iterator;
    typedef typename hypersurface_t::hypersurface_skeleton_t hypersurface_skeleton_t;
    
    typedef adjacent_vertex_pred<hypersurface_skeleton_t> adj_vert_pred_t;
    typedef boost::filtered_graph<hypersurface_skeleton_t,boost::keep_all,adj_vert_pred_t> filtered_graph_t;
    typedef boost::graph_traits<filtered_graph_t> filtered_traits;
    
    const hypersurface_skeleton_t& hyp_skel = hypersurface.hypersurface_skeleton;
    v_iterator v_iter, v_end;
    for(boost::tie(v_iter,v_end) = boost::vertices(hyp_skel); v_iter != v_end; ++v_iter) {
        //define subgraph containing all adjacent vertices -> link of *iter
        adj_vert_pred_t pred(*v_iter,hyp_skel);
        filtered_graph_t f_hyp_skel(hyp_skel,boost::keep_all(),pred);
        
        //count the vertices, should be >=4 (tetrahedron homeomorphic to 2-sphere)
        unsigned short verts = pred.num_adjacent_vertices();
        if(verts < 4)
            return false;
        
        //count edges in link of *iter
        typename filtered_traits::edge_iterator e_iter, e_end;
        unsigned short edges = 0;
        for(boost::tie(e_iter,e_end) = boost::edges(f_hyp_skel); e_iter != e_end; ++e_iter) {
            edges++;
        }
        if(edges < 6)
            return false;
        
        //count connected components in link of *iter, should be one
        std::vector<int> components(pred.num_adjacent_vertices()+1);
        if(boost::connected_components(f_hyp_skel,&components[0]) != 1)
            return false;
        
        //count triangles in link of *iter
        typename hypersurface_t::vector_set_t::const_iterator t_iter, t_end = hypersurface.triangle_set.end();
        unsigned short triangles = 0;
        for(t_iter = hypersurface.triangle_set.begin(); t_iter != t_end; ++t_iter) {
            typename hypersurface_t::vertex_container_t triangle_verts = *t_iter;
            bool isInLink = true;
            for(auto t_v_iter = triangle_verts.begin(); (t_v_iter != triangle_verts.end() && isInLink); ++t_v_iter) {
                isInLink = isInLink && pred(*t_v_iter);
            }
            if(isInLink)
                triangles++;
        }
        
        //compute euler characteristic
        return (verts - edges + triangles) == 2;
    } 
    
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::generateTriangulation() {
    // color the hypersurface to get independent vertex sets
    this->colorHypersurface();
    // generate 3 hypersurfaces / foliations, 0 -> 1 initial conditions, 1 -> 2 newton iteration
    for(hyper_index_t hyper_index = 0; hyper_index < 2; ++hyper_index) {
        // evolve each independent set of vertices one after another 
        // could be parallelized, but boost graph (vertex/edge insertions) is not thread safe (but there is a parallel graph lib)!
        for(vertex_num_t color = 0; color < colors; ++color) {
            vertex_color_pred<skeleton_t> pred(color,hyper_index,skeleton);
            typedef boost::filtered_graph<skeleton_t,boost::keep_all,vertex_color_pred<skeleton_t>> filtered_skeleton_t;
            filtered_skeleton_t filtered_skeleton(skeleton,boost::keep_all(),pred);
            typename boost::graph_traits<filtered_skeleton_t>::vertex_iterator iter, end;
            for(boost::tie(iter,end) = boost::vertices(filtered_skeleton); iter!=end; ++iter) {
                this->tentMove(*iter);
            }
        }
    }
    
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::tentMove(vertex_t vertex) {
    hyper_index_t current_index = skeleton[vertex].hyper_index;
    vertex_t next = boost::add_vertex(vertex_property_t(vertex,skeleton[vertex].color,current_index+1),skeleton);

    skeleton[vertex].next = next;
    skeleton[next].last = vertex;
    
    // collect "adjacent" vertices, to circumvent any iterator invalidation upon edge insertions
    std::vector<vertex_t> adj_vertices;
    adjacent_vertices_iterator iter, end;
    for(boost::tie(iter,end) = boost::adjacent_vertices(vertex,skeleton); iter != end; ++iter) {
        vertex_property_t& vert_prop = skeleton[*iter];
        if(vert_prop.hyper_index == current_index) {
            if(!vert_prop.evolved)
                adj_vertices.push_back(*iter);
            else {
                if(vert_prop.next == skeleton_descriptors::null_vertex())
                    throw std::logic_error("Tent move discoverd null vertex.");
                adj_vertices.push_back(vert_prop.next);
            }
        }
    }
    // connect vertices, create "tent"-like structure in underlying graph
    for(auto in_iter = adj_vertices.begin(); in_iter != adj_vertices.end(); ++in_iter) {
        boost::add_edge(*in_iter,next,edge_property_t(1.0),skeleton);
    }
    boost::add_edge(vertex,next,edge_property_t(-1.0),skeleton);
    skeleton[vertex].evolved = true;
    
    // start to identify simplices, create local filtered graph of new "tent"
    adj_vertices.push_back(vertex);
    vertex_vector_pred<skeleton_t> pred(adj_vertices);
    typedef boost::filtered_graph<skeleton_t,boost::keep_all,vertex_vector_pred<skeleton_t>> filtered_skeleton_t;
    filtered_skeleton_t filtered_skeleton(skeleton,boost::keep_all(),pred);
    
    // search for existing tetrahedra in new "tent", they form the bases for the new pentachora
    std::vector<vertex_container_t> tetrahedra;
    typename boost::graph_traits<filtered_skeleton_t>::edge_iterator e_iter, e_end;
    for(boost::tie(e_iter,e_end) = boost::edges(filtered_skeleton); e_iter != e_end; ++e_iter) {
        tetrahedron_container_t& e_tetras = filtered_skeleton[*e_iter].tetrahedrons;
        for(auto tetra_iter = e_tetras.begin(); tetra_iter != e_tetras.end(); ++tetra_iter) {
            tetrahedron_t& tetra = tetra_iter->get();                    
            bool isInTent = true;
            vertex_container_t temp_tetra;
            for(auto in_iter = tetra.edges.begin(); (in_iter != tetra.edges.end()) && isInTent; ++in_iter) {
                vertex_t src = boost::source(*in_iter,skeleton);
                vertex_t tgt = boost::target(*in_iter,skeleton);
                isInTent = (std::find(adj_vertices.begin(),adj_vertices.end(),src) != adj_vertices.end()) 
                        && (std::find(adj_vertices.begin(),adj_vertices.end(),tgt) != adj_vertices.end());
                if(isInTent) {
                    if(std::find(temp_tetra.begin(),temp_tetra.end(),src) == temp_tetra.end()) {
                        temp_tetra.push_back(src);
                    }
                    if(std::find(temp_tetra.begin(),temp_tetra.end(),tgt) == temp_tetra.end()) {
                        temp_tetra.push_back(tgt);
                    }
                }
            }
            if(isInTent) {
                tetrahedra.push_back(temp_tetra);                
            }
        }
        
    }
    // create pentachora on the tetrahedral bases
    std::vector<vertex_container_t> pentachora(std::move(tetrahedra));
    for(auto iter = pentachora.begin(); iter != pentachora.end(); ++ iter) {
        iter->push_back(next);
        std::sort(iter->begin(),iter->end());
        this->addSimplices(*iter);
    }
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::addSimplices(const vertex_container_t& vertices) {
    typedef edge_container_t triangle_faces_t;
    typedef std::pair<edge_container_t,std::vector<triangle_faces_t>> tetrahedron_faces_t;
    typedef std::pair<edge_container_t,std::vector<tetrahedron_faces_t>> pentachoron_faces_t;
    
    // create subsets of vertices
    pentachoron_faces_t subsets = edge_subsets_generator()(vertices,skeleton);
    
    for(auto iter = subsets.second.begin(); iter != subsets.second.end(); ++iter) {
        tetrahedron_t& tetrahedron = this->getOrCreateSimplex(tetrahedron_map,iter->first);
        for(auto tri_iter = iter->second.begin(); tri_iter != iter->second.end(); ++tri_iter) {
            triangle_t& triangle = this->getOrCreateSimplex(triangle_map,*tri_iter);
            this->addSimplexToContainer(triangle.tetrahedrons,tetrahedron);
            this->addSimplexToContainer(tetrahedron.triangles,triangle);
            
            for(auto in_iter = triangle.edges.begin(); in_iter != triangle.edges.end(); ++in_iter) {
                // auto e = boost::edge(boost::source(*in_iter,skeleton),boost::target(*in_iter,skeleton),skeleton).first;
                this->addSimplexToContainer(skeleton[*in_iter].triangles,triangle);
            }
            for(auto in_iter = tetrahedron.edges.begin(); in_iter != tetrahedron.edges.end(); ++in_iter) {
                // auto e = boost::edge(boost::source(*in_iter,skeleton),boost::target(*in_iter,skeleton),skeleton).first;
                this->addSimplexToContainer(skeleton[*in_iter].tetrahedrons,tetrahedron);
            }
            
            if(subsets.first.size() == 10) {
                pentachoron_t& pentachoron = this->getOrCreateSimplex(pentachoron_map,subsets.first);
                this->addSimplexToContainer(pentachoron.triangles,triangle);
                this->addSimplexToContainer(pentachoron.tetrahedrons,tetrahedron);
                // this->addSimplexToContainer(triangle.tetrahedrons,tetrahedron);
                this->addSimplexToContainer(triangle.pentachorons,pentachoron);
                // this->addSimplexToContainer(tetrahedron.triangles,triangle);
                this->addSimplexToContainer(tetrahedron.pentachorons,pentachoron);
                
                for(auto in_iter = pentachoron.edges.begin(); in_iter != pentachoron.edges.end(); ++in_iter) {
                    auto e = boost::edge(boost::source(*in_iter,skeleton),boost::target(*in_iter,skeleton),skeleton).first;
                    this->addSimplexToContainer(skeleton[e].pentachorons,pentachoron);
                }
            }
        }
    }
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::colorHypersurface() {    
    boost::vector_property_map<vertex_t> order_map;
    boost::smallest_last_vertex_ordering(skeleton,order_map);
    
    auto color_map(boost::get(&vertex_property_t::color,skeleton));
    this->colors = boost::sequential_vertex_coloring(skeleton,order_map,color_map);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::copySimplices(const hypersurface_t& hypersurface) {    
    const typename hypersurface_t::vector_set_t& hypsurf_tetras = hypersurface.tetrahedron_set;
    for(auto iter = hypsurf_tetras.begin(); iter != hypsurf_tetras.end(); ++iter) {
        this->addSimplices(*iter);        
    }
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
const typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_property_t& 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getEdgeProperty(const edge_t edge) const {
    return skeleton[edge];
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
const typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_property_t& 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getVertexProperty(const vertex_t vertex) const {
    return skeleton[vertex];
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setEdgeLengthSq(const edge_t edge, double length_sq) {
    skeleton[edge].length_sq = length_sq;            
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setVertexWeight(const vertex_t v, const vertex_weight_t w) {
    skeleton[v].weight = w;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setEdgeWeight(const edge_t e, const edge_weight_t w) {
    skeleton[e].weight = w;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setTriangleWeight(const triangle_t t, const triangle_weight_t w) {
    triangle_map.at(t).second.weight = w;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setTetrahedronWeight(const tetrahedron_t t, const tetrahedron_weight_t w) {
    tetrahedron_map.at(t).second.weight = w;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setPentachoronWeight(const pentachoron_t p, const pentachoron_weight_t w) {
    pentachoron_map.at(p).second.weight = w;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_t 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getSource(const edge_t& e) const {
    return std::min(boost::source(e,skeleton),boost::target(e,skeleton));
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_t 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getTarget(const edge_t& e) const {
    return std::max(boost::source(e,skeleton),boost::target(e,skeleton));
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_iterator, typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_iterator> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getEdges() const {
    return boost::edges(skeleton);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_iterator, typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_iterator> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getVertices() const {
    return boost::vertices(skeleton);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::const_triangle_iterator, typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::const_triangle_iterator> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getTriangles() const {
    return std::make_pair(triangle_map.begin(),triangle_map.end());
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::const_tetrahedron_iterator, typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::const_tetrahedron_iterator> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getTetrahedra() const {
    return std::make_pair(tetrahedron_map.begin(),tetrahedron_map.end());
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::const_pentachoron_iterator, typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::const_pentachoron_iterator> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getPentachora() const {
    return std::make_pair(pentachoron_map.begin(),pentachoron_map.end());
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_t, bool> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getEdge(vertex_t v, vertex_t u) const {
    return boost::edge(v,u,skeleton);            
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::adjacent_vertices_iterator, typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::adjacent_vertices_iterator> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getAdjacentVertices(const vertex_t vertex) const {
    return boost::adjacent_vertices(vertex,skeleton);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_num_t
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getNumColors() const {
    return colors;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
unsigned long triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getNumVertices() const {
    return boost::num_vertices(skeleton);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
unsigned long triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getNumEdges() const {
    return boost::num_edges(skeleton);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
unsigned long triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getNumTetrahedra() const {
    return tetrahedron_map.size();
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
unsigned long triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getNumTriangles() const {
    return triangle_map.size();
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
unsigned long triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getNumPentachora() const {
    return pentachoron_map.size();
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
std::pair<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::hyper_index_t,typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_num_t> 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getInputIndices(const vertex_t vertex) const {
    hyper_index_t h = this->getVertexProperty(vertex).hyper_index;
    vertex_t v0 = this->getVertexAt(0,vertex);
    return std::make_pair(h,vertex_index_bimap.right.find(v0)->second);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_t 
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getVertex(hyper_index_t h, vertex_num_t index) const {
    vertex_t v0 = vertex_index_bimap.left.find(index)->second;
    return this->getVertexAt(h,v0);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::setEvolved(const vertex_t vertex) {
    skeleton[vertex].evolved = true;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
template<typename vertex_predicate_t>
const boost::filtered_graph<typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::skeleton_t,boost::keep_all,vertex_predicate_t>
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getFilteredSkeleton(const vertex_predicate_t& vertex_predicate) {
    return boost::filtered_graph<skeleton_t,boost::keep_all,vertex_predicate_t>(skeleton,boost::keep_all(),vertex_predicate);
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
const typename triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::skeleton_t&
triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::getSkeleton() {
    return skeleton;
}

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
void triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::swapHypersurfaces() {
    for(auto iter = this->getEdges().first; iter != this->getEdges().second; ++iter) {
        vertex_t src = this->getSource(*iter);
        vertex_t tgt = this->getTarget(*iter);
        if((skeleton[src].hyper_index > 0) && (skeleton[tgt].hyper_index > 0)) {
            if(!this->getEdge(skeleton[src].last,skeleton[tgt].last).second)
                throw std::logic_error("Trying to swap to non-existent edge after Sorkin step.");
            edge_t e_lower = this->getEdge(skeleton[src].last,skeleton[tgt].last).first;
            skeleton[e_lower] = skeleton[*iter];
        }
    }
    for(auto iter = this->getVertices().first; iter != this->getVertices().second; ++iter) {
        vertex_property_t& v_prop = skeleton[*iter];
        if(v_prop.hyper_index > 0) {
            skeleton[v_prop.last] = v_prop; 
            if(v_prop.hyper_index == 2)
                v_prop.evolved = false;
        }
    }
    simplex_swap_pred<skeleton_t> swap(this->skeleton);
    for(auto iter = triangle_map.begin(); iter != triangle_map.end(); ++iter) {
        iter->second.edges;
        if(std::find_if(iter->second.edges.begin(),iter->second.edges.end(),swap) == iter->second.edges.end()) {
            edge_container_t swap_simplex_edges;
            for(auto swap_iter = iter->second.edges.begin(); swap_iter != iter->second.edges.end(); ++swap_iter) {
                vertex_t last_src = skeleton[this->getSource(*swap_iter)].last;
                vertex_t last_tgt = skeleton[this->getTarget(*swap_iter)].last;
                
                auto edge_pair = this->getEdge(last_src,last_tgt);
                if(!edge_pair.second)
                    throw std::logic_error("Trying to swap to non-existent triangle after Sorkin step.");
                swap_simplex_edges.push_back(edge_pair.first);
            }
            triangle_map.find(triangle_t(swap_simplex_edges))->second = iter->second;
        }
    }
    for(auto iter = tetrahedron_map.begin(); iter != tetrahedron_map.end(); ++iter) {
        iter->second.edges;
        if(std::find_if(iter->second.edges.begin(),iter->second.edges.end(),swap) == iter->second.edges.end()) {
            edge_container_t swap_simplex_edges;
            for(auto swap_iter = iter->second.edges.begin(); swap_iter != iter->second.edges.end(); ++swap_iter) {
                vertex_t last_src = skeleton[this->getSource(*swap_iter)].last;
                vertex_t last_tgt = skeleton[this->getTarget(*swap_iter)].last;
                
                auto edge_pair = this->getEdge(last_src,last_tgt);
                if(!edge_pair.second)
                    throw std::logic_error("Trying to swap to non-existent tetrahedron after Sorkin step.");
                swap_simplex_edges.push_back(edge_pair.first);
            }
            tetrahedron_map.find(tetrahedron_t(swap_simplex_edges))->second = iter->second;
        }
    }
    for(auto iter = pentachoron_map.begin(); iter != pentachoron_map.end(); ++iter) {
        iter->second.edges;
        if(std::find_if(iter->second.edges.begin(),iter->second.edges.end(),swap) == iter->second.edges.end()) {
            edge_container_t swap_simplex_edges;
            for(auto swap_iter = iter->second.edges.begin(); swap_iter != iter->second.edges.end(); ++swap_iter) {
                vertex_t last_src = skeleton[this->getSource(*swap_iter)].last;
                vertex_t last_tgt = skeleton[this->getTarget(*swap_iter)].last;
                
                auto edge_pair = this->getEdge(last_src,last_tgt);
                if(!edge_pair.second)
                    throw std::logic_error("Trying to swap to non-existent pentachoron after Sorkin step.");
                swap_simplex_edges.push_back(edge_pair.first);
            }
            pentachoron_map.find(pentachoron_t(swap_simplex_edges))->second = iter->second;
        }
    }
}

template<>
void triangulator_t<empty_weight, empty_weight, empty_weight, empty_weight, empty_weight>::swapHypersurfaces() {
    // copy 1st hypersurface first (h = 1)
    vertex_hyper_index_pred<skeleton_t> pred(1,skeleton);
    auto filtered_skeleton = this->getFilteredSkeleton(pred);
    for(auto iter = boost::edges(filtered_skeleton).first; iter != boost::edges(filtered_skeleton).second; ++iter) { 
        vertex_t src = this->getSource(*iter);
        vertex_t tgt = this->getTarget(*iter);
        edge_t e_lower = this->getEdge(skeleton[src].last,skeleton[tgt].last).first;
        skeleton[e_lower] = skeleton[*iter];
    }
    
    // copy rest
    for(auto iter = this->getEdges().first; iter != this->getEdges().second; ++iter) {
        vertex_t src = this->getSource(*iter);
        vertex_t tgt = this->getTarget(*iter);
        
        if(skeleton[tgt].hyper_index == 2) {
            edge_t e_lower = this->getEdge(skeleton[src].last,skeleton[tgt].last).first;
            skeleton[e_lower] = skeleton[*iter];
        }
    }
    
    // same for vertices
    
    for(auto iter = boost::vertices(filtered_skeleton).first; iter != boost::vertices(filtered_skeleton).second; ++iter) { 
        vertex_property_t& v_prop = skeleton[*iter];
        skeleton[v_prop.last] = v_prop;
    }
    
    for(auto iter = this->getVertices().first; iter != this->getVertices().second; ++iter) {
        vertex_property_t& v_prop = skeleton[*iter];
        if(v_prop.hyper_index == 2) {
            skeleton[v_prop.last] = v_prop; 
            v_prop.evolved = false;
        }
    }
}

#endif	/* TRIANGULATOR_IMPL_H */


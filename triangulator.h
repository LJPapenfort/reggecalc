/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef TRIANGULATOR_H
#define	TRIANGULATOR_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/bimap.hpp>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "regge_generator.h"

class empty_weight {};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t {
    public:
        // initial hypersurface triangulation class
        class hypersurface_t {
            friend class triangulator_t<vertex_weight_t,edge_weight_t,triangle_weight_t,tetrahedron_weight_t,pentachoron_weight_t>;
            public:
                typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> hypersurface_skeleton_t;
                typedef boost::graph_traits<hypersurface_skeleton_t> traits;
                typedef typename traits::vertex_descriptor vertex_t;
                typedef typename traits::vertices_size_type num_vertex_size_t;
                typedef typename std::vector<vertex_t> vertex_container_t;
                typedef typename std::vector<vertex_container_t> tetrahedron_container_t;
                typedef std::unordered_set<vertex_container_t,vector_hash<vertex_t>> vector_set_t;
            private:
                typedef boost::bimap<num_vertex_size_t,num_vertex_size_t> vertex_index_bimap_t;
                typedef vertex_index_bimap_t::value_type vertex_indices_t;
                hypersurface_skeleton_t hypersurface_skeleton;
                vector_set_t tetrahedron_set;
                vector_set_t triangle_set;
                vertex_index_bimap_t vertex_index_bimap;
                void clear();
            public:
                hypersurface_t() { }
                hypersurface_t(tetrahedron_container_t& vec);
                void addTetrahedron(const vertex_t v1, const vertex_t v2, const vertex_t v3, const vertex_t v4);
                void addTetrahedron(vertex_container_t& vec);
                void printSimplices();
        };
        
        // data structures
        class simplex_t;
        class vertex_property_t;
        class edge_property_t;
        class triangle_t;
        class tetrahedron_t;
        class pentachoron_t;
        
        // functor classes
        class simplex_t_hash;
        class edge_subsets_generator;
        class edge_t_equality;
        class edge_t_comparison;
        
        // geodesic base class, to deliver geodesic lengths / line elements
        
        class geodesics_t;
        
    public:     
        // 1-skeleton ^= underlying graph (vertices and edges) of the triangulation types
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, vertex_property_t, edge_property_t> skeleton_t;
        typedef boost::graph_traits<skeleton_t> skeleton_descriptors;
        typedef typename skeleton_descriptors::vertex_descriptor vertex_t;
        typedef typename skeleton_descriptors::edge_descriptor edge_t;
        typedef typename skeleton_descriptors::adjacency_iterator adjacent_vertices_iterator;
        typedef typename skeleton_descriptors::edge_iterator edge_iterator;
        typedef typename skeleton_descriptors::vertex_iterator vertex_iterator;
        
        // local vertex / edge property container types
        typedef std::vector<edge_t> edge_container_t;
        typedef std::vector<vertex_t> vertex_container_t;
        typedef std::vector<std::reference_wrapper<triangle_t>> triangle_container_t;
        typedef std::vector<std::reference_wrapper<tetrahedron_t>> tetrahedron_container_t;
        typedef std::vector<std::reference_wrapper<pentachoron_t>> pentachoron_container_t;
        
        // global simplex map types
        typedef std::unordered_map<triangle_t,triangle_t,simplex_t_hash> triangle_map_t;
        typedef std::unordered_map<tetrahedron_t,tetrahedron_t,simplex_t_hash> tetrahedron_map_t;
        typedef std::unordered_map<pentachoron_t,pentachoron_t,simplex_t_hash> pentachoron_map_t;
        typedef typename triangle_map_t::const_iterator const_triangle_iterator;
        typedef typename tetrahedron_map_t::const_iterator const_tetrahedron_iterator;
        typedef typename pentachoron_map_t::const_iterator const_pentachoron_iterator;
        
        // vertex index bimap
        typedef typename skeleton_descriptors::vertices_size_type vertex_num_t;
        typedef boost::bimap<vertex_num_t,vertex_num_t> vertex_index_bimap_t;
        typedef unsigned short hyper_index_t;
        
        
        triangulator_t(hypersurface_t& hypersurface, const geodesics_t& geodesics);
        ~triangulator_t(); 
        
        void swapHypersurfaces();
                
        void setEdgeLengthSq(const edge_t edge, double length_sq);
        void setEvolved(const vertex_t vertex);
        void setVertexWeight(const vertex_t v, const vertex_weight_t w);
        void setEdgeWeight(const edge_t e, const edge_weight_t w);
        void setTriangleWeight(const triangle_t t, const triangle_weight_t w);
        void setTetrahedronWeight(const tetrahedron_t t, const tetrahedron_weight_t w);
        void setPentachoronWeight(const pentachoron_t p, const pentachoron_weight_t w);
        
        template<typename vertex_predicate_t> const boost::filtered_graph<skeleton_t,boost::keep_all,vertex_predicate_t> getFilteredSkeleton(const vertex_predicate_t& vertex_predicate);
        const skeleton_t& getSkeleton();
        
        const edge_property_t& getEdgeProperty(const edge_t edge) const;
        const vertex_property_t& getVertexProperty(const vertex_t vertex) const;
        vertex_t getSource(const edge_t& e) const;
        vertex_t getTarget(const edge_t& e) const;
        std::pair<edge_iterator,edge_iterator> getEdges() const;
        std::pair<edge_t, bool> getEdge(vertex_t v, vertex_t u) const;
        std::pair<vertex_iterator,vertex_iterator> getVertices() const;
        std::pair<const_triangle_iterator,const_triangle_iterator> getTriangles() const;
        std::pair<const_tetrahedron_iterator,const_tetrahedron_iterator> getTetrahedra() const;
        std::pair<const_pentachoron_iterator,const_pentachoron_iterator> getPentachora() const;
        std::pair<adjacent_vertices_iterator,adjacent_vertices_iterator> getAdjacentVertices(const vertex_t vertex) const;
        unsigned long getNumVertices() const;
        unsigned long getNumEdges() const;
        unsigned long getNumTetrahedra() const;
        unsigned long getNumTriangles() const;
        unsigned long getNumPentachora() const;
        vertex_num_t getNumColors() const;
        std::pair<hyper_index_t,vertex_num_t> getInputIndices(const vertex_t vertex) const;
        vertex_t getVertex(hyper_index_t h, vertex_num_t index) const;
        template<typename map_t> static typename map_t::mapped_type& getOrCreateSimplex(map_t& map, const edge_container_t& edges);
        template<typename type_t> static void addSimplexToContainer(std::vector<std::reference_wrapper<type_t>>& container, type_t& type);
    private:
        skeleton_t skeleton;
        triangle_map_t triangle_map;
        tetrahedron_map_t tetrahedron_map;
        pentachoron_map_t pentachoron_map;
        vertex_index_bimap_t vertex_index_bimap;
        
        
        vertex_num_t colors;
        vertex_t null_vertex;
        bool checkHypersurface(const hypersurface_t& hypersurface) const;
        void copySimplices(const hypersurface_t& hypersurface);
        void colorHypersurface();
        void generateTriangulation();
        void tentMove(vertex_t vertex);
        void addSimplices(const vertex_container_t& vertices);
        void setInitialLengths(const geodesics_t& geodesics);
        vertex_t getVertexAt(hyper_index_t h, vertex_t v) const;
};

#include "hypersurface.h"
#include "triangulator_functors_impl.h"
#include "triangulator_types_impl.h"
#include "triangulator_impl.h"
#endif	/* TRIANGULATOR_H */


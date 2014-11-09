/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef TRIANGULATOR_TYPES_IMPL_H
#define	TRIANGULATOR_TYPES_IMPL_H
#include "triangulator.h"
//#include <random>
//#include <chrono>


template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::geodesics_t {
    public:
        virtual double operator()(const hyper_index_t h1, const vertex_num_t index1, const hyper_index_t h2, const vertex_num_t index2) const = 0;
};


template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::simplex_t {
    protected:
        virtual int num_edges() const  {
            return 0; // by definition
        } 
    public:    
        edge_container_t edges;
        simplex_t(const edge_container_t& input_edges) {
            typedef typename edge_container_t::const_iterator c_iterator;
            
            int num_edges = this->num_edges();
            
            /* NOT WORKING
            if(!(input_edges.size() == num_edges)) {
                std::cout << num_edges;
                throw std::logic_error("Simplex initialized with false number of edges!");
            }
             */
            //edges.reserve(num_edges);
            for(c_iterator iter = input_edges.begin(); iter != input_edges.end(); ++iter) {
                edges.push_back(*iter);
            }
            std::sort(edges.begin(),edges.end(),edge_t_comparison());
        }
        friend bool operator==(const simplex_t& simplex_lhs,const simplex_t& simplex_rhs) {
            return simplex_lhs.edges == simplex_rhs.edges;
        }
        friend std::ostream& operator << (std::ostream &os, simplex_t& simplex) {
            os << "[";
            for(auto iter = simplex.edges.begin(); iter != simplex.edges.end(); ++iter) {
                os << *iter;
                if(std::next(iter) != simplex.edges.end())
                    os << ",";
            }
            os << "]";
            return os;
        }
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::vertex_property_t {
    public:
        vertex_t next;
        vertex_t last;
        bool evolved;
        vertex_num_t color;
        hyper_index_t hyper_index;
        vertex_weight_t weight;
        vertex_property_t() : vertex_property_t(skeleton_descriptors::null_vertex(),0,0) {}
        vertex_property_t(vertex_t last, vertex_num_t color, hyper_index_t hyper_index) : 
                next(skeleton_descriptors::null_vertex()) , evolved(false) , last(last) , color(color) , hyper_index(hyper_index)  {}
        vertex_property_t& operator=(const boost::no_property& no_prop) {
            next = skeleton_descriptors::null_vertex();
            evolved = false;
            last = next;
            color = 0;
            hyper_index = 0;
            return *this;
        } 
        vertex_property_t& operator=(const vertex_property_t& vertex_property) {
            weight = vertex_property.weight;
        }
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::edge_property_t {
    public:
        triangle_container_t triangles;
        tetrahedron_container_t tetrahedrons;
        pentachoron_container_t pentachorons;
        double length_sq;
        edge_weight_t weight;
        edge_property_t() : length_sq(0) {}
        edge_property_t(double length_sq) : length_sq(length_sq) {}
        edge_property_t& operator=(const boost::no_property& no_prop) {
            /*
            // random distributed length for testing purposes
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            std::uniform_real_distribution<double> distribution(0.0,1.0);
            length_sq = distribution(generator);
             */
            length_sq = 1.0d;
            return *this;
        }
        
        edge_property_t& operator=(const edge_property_t& e_prop) {
            length_sq = e_prop.length_sq;
            weight = e_prop.weight;
            return *this;
        } 
        edge_property_t(const edge_property_t& e_prop) : triangles(e_prop.triangles), tetrahedrons(e_prop.tetrahedrons), 
                pentachorons(e_prop.pentachorons), length_sq(e_prop.length_sq), weight(e_prop.weight) {}
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::triangle_t
        : public simplex_t {    
    private:
        int num_edges() const {
            return 3;
        }
    public:
        tetrahedron_container_t tetrahedrons;
        pentachoron_container_t pentachorons;
        triangle_weight_t weight;
        triangle_t& operator=(const triangle_t& triangle) {
            weight = triangle.weight;
        }
        using simplex_t::simplex_t;
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::tetrahedron_t
        : public simplex_t {
    private:
        int num_edges() const {
            return 6;
        }
    public:
        triangle_container_t triangles;
        pentachoron_container_t pentachorons;
        tetrahedron_weight_t weight;
        tetrahedron_t& operator=(const tetrahedron_t& tetrahedron) {
            weight = tetrahedron.weight;
        }
        using simplex_t::simplex_t;
};

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t>::pentachoron_t
        : public simplex_t {
    private:
        int num_edges() const {
            return 10;
        }
    public:
        triangle_container_t triangles;
        tetrahedron_container_t tetrahedrons;
        pentachoron_weight_t weight;
        pentachoron_t& operator=(const pentachoron_t& pentachoron) {
            weight = pentachoron.weight;
        }
        using simplex_t::simplex_t;
};

#endif	/* TRIANGULATOR_TYPES_IMPL_H */


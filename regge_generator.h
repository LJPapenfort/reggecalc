/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */


#ifndef REGGE_GENERATOR_H
#define	REGGE_GENERATOR_H

#include "triangulator.h"
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>

class empty_weight;

template<typename vertex_weight_t = empty_weight, typename edge_weight_t = empty_weight, typename triangle_weight_t = empty_weight, 
        typename tetrahedron_weight_t = empty_weight, typename pentachoron_weight_t = empty_weight>
class triangulator_t;

template<typename vertex_weight_t = empty_weight, typename edge_weight_t = empty_weight, typename triangle_weight_t = empty_weight, 
        typename tetrahedron_weight_t = empty_weight, typename pentachoron_weight_t = empty_weight>
class regge_functor_t {
    private:
        int getEdgeIndex(const int i, const int j) const {
            // helper to convert matrix indices to vector index
            return j+4*i-i*(i+1)/2-1;
        }
        inline double kahanSummation(double& compensation, const double summand, const double sumd) const {
            double corrected = summand - compensation;
            double temp_sum = sumd + corrected;
            compensation = (temp_sum - sumd) - corrected;
            return temp_sum;
        }
    protected:
        typedef triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> _triangulator_t;
        typedef typename _triangulator_t::simplex_t simplex_t;
        typedef typename _triangulator_t::pentachoron_t pentachoron_t;
        typedef typename _triangulator_t::edge_t edge_t;
        typedef typename _triangulator_t::edge_property_t edge_property_t;
        typedef typename _triangulator_t::edge_container_t edge_container_t;
        typedef typename _triangulator_t::vertex_t vertex_t;
        typedef typename _triangulator_t::vertex_property_t vertex_property_t;
        typedef typename _triangulator_t::vertex_container_t vertex_container_t;
        
        typedef Eigen::Matrix4d Matrix4d;
        
        _triangulator_t& triangulator;
        
        virtual inline double sumUpKahan(const std::vector<double>& doubles) const {
            double result = 0;
            double compensation = 0;
            for(auto iter = doubles.begin(); iter != doubles.end(); ++iter) {
                // Kahan summation, O(const) truncation error bound, but more operations
                result = this->kahanSummation(compensation,*iter,result); 
            }
            return result;
        }  
        virtual inline double sumUpNaive(const std::vector<double>& doubles) const {
            double result = 0;
            for(auto iter = doubles.begin(); iter != doubles.end(); ++iter) {
                // simple, but maybe numerically unstable (at least for lots of summands)
                result += *iter;
            }
            return result;
        }   
        virtual inline double sumUp(const std::vector<double>& doubles) const {
            return this->sumUpNaive(doubles);
            //return this->sumUpKahan(doubles);
        }
        virtual Matrix4d getMetric(const edge_container_t& edges) const {
            // calculate local metric from edge lenghts
            Matrix4d metric;
            for(int i = 1; i<=4 ; ++i) {
                for(int j = 1; j<=4; j++) {
                    double len_sq1 = triangulator.getEdgeProperty(edges.at(this->getEdgeIndex(0,i))).length_sq;
                    double len_sq2 = triangulator.getEdgeProperty(edges.at(this->getEdgeIndex(0,j))).length_sq;
                    double len_sq3 = 0;
                    if(i != j) 
                        len_sq3 = triangulator.getEdgeProperty(edges.at(this->getEdgeIndex(std::min(i,j),std::max(i,j)))).length_sq;
                    metric(i-1,j-1) = 1.0d/2.0d*(len_sq1+len_sq2-len_sq3); 
                }
            }
            return metric;
        }
        virtual void transformMetric(Matrix4d& metric, const edge_container_t& edges, const simplex_t& to_origin) const {
            // collect local metric vertices, e.g. (0,1,2,3,4)
            vertex_container_t metric_vertices;
            vertex_t origin = triangulator.getSource(*edges.begin());
            
            metric_vertices.push_back(origin);
            
            for(auto iter = edges.begin(); triangulator.getSource(*iter) == origin; ++iter) {
                metric_vertices.push_back(triangulator.getTarget(*iter));
            }
            
            // collect simplex vertices, e.g. (x,y,z)
            vertex_t new_origin = triangulator.getSource(*to_origin.edges.begin());
            vertex_container_t simplex_vertices;
            simplex_vertices.push_back(new_origin);
            for(auto iter = to_origin.edges.begin(); triangulator.getSource(*iter) == new_origin; ++iter) {
                simplex_vertices.push_back(triangulator.getTarget(*iter));
            }
            
            // look for positions of vertices in the standard frame
            std::vector<int> vertex_positions;
            vertex_positions.reserve(metric_vertices.size());
            for(auto iter = simplex_vertices.begin(); iter != simplex_vertices.end(); ++iter) {
                auto find_it = std::find(metric_vertices.begin(),metric_vertices.end(),*iter);
                if(find_it == metric_vertices.end())
                    throw std::logic_error("Could not locate position of vertex in standard frame.");
                // std::distance is O(1) for random access containers
                vertex_positions.push_back(std::distance(metric_vertices.begin(),find_it));
            }
            
            // sort out identity transformations
            if(vertex_positions.back() == vertex_positions.size()-1) {
                return ;
            }
            
            // add missing vertex positions
            int temp = -1;
            for(auto iter = vertex_positions.begin(); (iter != vertex_positions.end()) ; ++iter) {
                while(*iter > ++temp) {
                    vertex_positions.push_back(temp);
                }      
            }
            while(vertex_positions.size() < metric_vertices.size()) {
                vertex_positions.push_back(vertex_positions.size());
            }
            
            // generate transformation matrix from vertex positions
            Matrix4d trans_matrix;
            for(int i = 0; i<4; ++i) {
                for(int k = 0; k<4; ++k) {
                    if(i+1 == vertex_positions.front())
                        trans_matrix(i,k) = -1.0d;
                    else
                        trans_matrix(i,k) = 0.0d;
                }
            }
            int j = 0;
            for(auto iter = ++vertex_positions.begin(); iter != vertex_positions.end(); ++iter) {
                if(*iter != 0)
                    trans_matrix(*iter-1,j) = 1.0d;
                ++j;
            }
            
            // transform metric
            metric = trans_matrix.transpose() * metric * trans_matrix;
        }
        inline double PI() const {
            return boost::math::constants::pi<double>();
        }
    public:
        virtual double operator()(const edge_t& edge) const = 0;
        virtual double partial_derivative(const edge_t& edge, const edge_t& derivative_edge) {
             // numerical differentiation, central difference quotient as default
            // should be overloaded!
            
            //std::cout << "Computing derivative of " << edge << " by " << derivative_edge << std::endl;
            double central_len_sq = this->triangulator.getEdgeProperty(derivative_edge).length_sq;
            double h = std::sqrt(std::numeric_limits<double>::epsilon())*central_len_sq;
            this->triangulator.setEdgeLengthSq(derivative_edge,central_len_sq+h);
            double temp = this->operator ()(edge);
            this->triangulator.setEdgeLengthSq(derivative_edge,central_len_sq-h);
            temp = temp - this->operator ()(edge);
            temp = temp/(2*h);
            
            this->triangulator.setEdgeLengthSq(derivative_edge,central_len_sq);

            return temp;
        }
        regge_functor_t(_triangulator_t& t) : triangulator(t) {};
};

template<typename vertex_weight_t = empty_weight, typename edge_weight_t = empty_weight, typename triangle_weight_t = empty_weight, 
        typename tetrahedron_weight_t = empty_weight, typename pentachoron_weight_t = empty_weight>
class regge_generator_t : public regge_functor_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> {
    private:
        typedef regge_functor_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> _regge_functor_t;
        typedef typename _regge_functor_t::edge_property_t edge_property_t;
        typedef typename _regge_functor_t::_triangulator_t _triangulator_t;
        typedef typename _regge_functor_t::edge_t edge_t;
        std::vector<_regge_functor_t*> functors;
    public:
        template<typename... regge_functors> regge_generator_t(_triangulator_t& t, regge_functors&... funcs) : _regge_functor_t(t) , functors{&funcs...} {}
        ~regge_generator_t() {};
        virtual double operator()(const edge_t& edge) const {
            std::vector<double> results;
            for(auto iter = functors.begin(); iter != functors.end(); ++iter) {
                double result = (*(*iter))(edge);
                results.push_back(result); 
            }
            return this->sumUp(results);
        }
        virtual double partial_derivative(const edge_t& edge, const edge_t& derivative_edge) {
            std::vector<double> results;
            for(auto iter = functors.begin(); iter != functors.end(); ++iter) {
                double result = (*iter)->partial_derivative(edge,derivative_edge);
                results.push_back(result); 
            }
            return this->sumUp(results);
        }
};

#include "regge_functor_special.h"
#endif	/* REGGE_GENERATOR_H */


/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef SORKIN_STEPPER_H
#define	SORKIN_STEPPER_H

#include <thread>
#include "triangulator.h"
#include "utilities.h"
#include "regge_generator.h"

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class sorkin_stepper {
        typedef triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> _triangulator_t;
        typedef typename _triangulator_t::skeleton_t skeleton_t;
        typedef typename _triangulator_t::vertex_t vertex_t;
        typedef typename _triangulator_t::vertex_property_t vertex_property_t;
        typedef typename _triangulator_t::edge_t edge_t;
        typedef typename _triangulator_t::edge_container_t edge_container_t;
        typedef typename _triangulator_t::vertex_num_t vertex_num_t;
        typedef regge_generator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> _regge_generator_t;
        
        typedef boost::filtered_graph<skeleton_t,boost::keep_all,vertex_hyper_index_pred<skeleton_t>> hyper_index_filtered_skeleton_t;
        
        _triangulator_t& triangulator;
        _regge_generator_t& regge_generator;

        double courant;
        unsigned threads;
        
        template<typename _shift_conditions_t> void computeRHS(const edge_container_t& equations, _shift_conditions_t& shift_conditions, Eigen::VectorXd& rhs) const {
            for(int i = 0; i<equations.size(); ++i) {
                rhs(i) = -regge_generator(equations.at(i));
            }
            const std::vector<double>& shift_rhs = shift_conditions();
            for(int i = 0; i<shift_rhs.size(); ++i) {
                rhs(i+equations.size()) = -shift_rhs.at(i);
            }
        }
        
        void updateEdgeLengthsSq(const edge_container_t& unknowns, const Eigen::VectorXd& dlength_sqs) {
            for(int i = 0; i<unknowns.size(); ++i) {
                triangulator.setEdgeLengthSq(unknowns.at(i),triangulator.getEdgeProperty(unknowns.at(i)).length_sq+dlength_sqs(i));
            }
        }
        
        template<typename _shift_conditions_t> void solve(vertex_t vertex, const hyper_index_filtered_skeleton_t& hyper_index_filtered_skeleton, double dt_sq, _shift_conditions_t shift_conditions) {
            edge_container_t unknowns;
            edge_container_t equations;
            
            auto adj_pair = boost::adjacent_vertices(vertex,hyper_index_filtered_skeleton);
            for(auto adj_iter = adj_pair.first; adj_iter != adj_pair.second; ++adj_iter) {
                const vertex_property_t adj_prop = triangulator.getVertexProperty(*adj_iter);
                if(adj_prop.evolved) {
                    auto unknown_pair = triangulator.getEdge(vertex,*adj_iter);
                    auto equation_pair = triangulator.getEdge(triangulator.getVertexProperty(vertex).last,*adj_iter);
                    if((!unknown_pair.second) || (!equation_pair.second))
                        throw std::logic_error("Encouted non-existing spacelike edge in sorkin-stepper->solve routine.");
                    unknowns.push_back(unknown_pair.first);
                    equations.push_back(equation_pair.first);
                }
                else {
                    auto unknown_pair = triangulator.getEdge(vertex,triangulator.getVertexProperty(*adj_iter).last);
                    auto equation_pair = triangulator.getEdge(triangulator.getVertexProperty(vertex).last,triangulator.getVertexProperty(*adj_iter).last);
                    if((!unknown_pair.second) || (!equation_pair.second))
                        throw std::logic_error("Encouted non-existing spacelike edge in sorkin-stepper->solve routine.");
                    unknowns.push_back(unknown_pair.first);
                    equations.push_back(equation_pair.first);                           
                }
            }
            // get / set lapse edge
            auto lapse_pair = triangulator.getEdge(triangulator.getVertexProperty(vertex).last,vertex);
            if(!lapse_pair.second)
                throw std::logic_error("Encouted non-existing timelike edge in sorkin-stepper->solve routine.");
            triangulator.setEdgeLengthSq(lapse_pair.first,-dt_sq);
            
            equations.push_back(lapse_pair.first);
            
            // set / get shift
            shift_conditions.setUnknowns(unknowns);
            /* Activate this if you want to do a LU scheme, omitting the Regge equations associated with the shift edges
            for(auto iter = shift_conditions.getConstraintEdges().begin(); iter != shift_conditions.getConstraintEdges().end(); ++iter) {
                bool found = false;
                for(int i = 0; !found && (i<unknowns.size()); ++i) {
                    if(unknowns.at(i) == *iter) {
                        found = true;
                        equations.erase(equations.begin()+i);
                    }
                }
            }
             */
             
            // init matrix and rhs of Ax = b
            int n = equations.size()+shift_conditions.getNumConditions();
            int m = unknowns.size();
            //std::cout << "Got n=" << n << " equations and m=" << m << " unknowns" << std::endl;
            Eigen::MatrixXd A(n,m);
            Eigen::VectorXd b(n);
            Eigen::VectorXd solution(m);
            // init solver
            Eigen::FullPivHouseholderQR<Eigen::MatrixXd> solver(n,m);
            
            // different possibilities, depending on the system of equations
            //Eigen::FullPivLU<Eigen::MatrixXd> solver(n,m);
            //Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(n,m);
            //Eigen::JacobiSVD<Eigen::MatrixXd> solver(n,m,Eigen::ComputeThinU | Eigen::ComputeThinV);

            // init rhs
            this->computeRHS(equations,shift_conditions,b);
            
            /*
            std::cout << "Got " << n << " equations and " << m << " unknowns." << std::endl;
            std::cout << "Initial RHS: ";
            for(int i = 0; i<n; ++i) {
                std::cout << b(i) << " ";
            }
            std::cout << std::endl;
             */
            
            int j = 1;
            //while(std::abs(b_norm - b.norm()) > 1.0e-4) {
            //while(!b.isZero(thres_i) || (j<10)) {
            while(j<10) {
                for(int i = 0; i < equations.size(); ++i) {
                    for(int j = 0; j < m; ++j) {
                        A(i,j) = regge_generator.partial_derivative(equations.at(i),unknowns.at(j));
                    }
                }
                for(int i = 0; i<m; ++i) {
                    const std::vector<double>& shift_partials = shift_conditions.partial_derivative(unknowns.at(i));
                    for(int j = 0; j<shift_partials.size(); ++j) {
                        A(j+equations.size(),i) = shift_partials.at(j);
                    }
                }
                if(A.hasNaN()) {
                    std::cout << A << std::endl;
                    throw std::logic_error("A has NaNs");
                }
                solver.compute(A);
                solution = solver.solve(b);
                if(solution.hasNaN())
                    throw std::logic_error("Solution has NaNs");
                //b_norm = b.norm();
                
                this->updateEdgeLengthsSq(unknowns,solution);
                this->computeRHS(equations,shift_conditions,b);
                if(b.hasNaN())
                    throw std::logic_error("RHS has NaNs");
                
                j++;
            }
            triangulator.setEvolved(vertex);
        }
        double getMinLength() {
            vertex_hyper_index_pred<skeleton_t> h_pred(1,triangulator.getSkeleton());
            auto filtered_skeleton = triangulator.getFilteredSkeleton(h_pred);
            double min_len_sq = triangulator.getEdgeProperty(*boost::edges(filtered_skeleton).first).length_sq;
            for(auto iter = boost::edges(filtered_skeleton).first; iter != boost::edges(filtered_skeleton).second; ++iter) {
                if(triangulator.getEdgeProperty(*iter).length_sq < min_len_sq)
                    min_len_sq = triangulator.getEdgeProperty(*iter).length_sq;
            }
            return min_len_sq;
        }
    public:
        // base class for shift conditions
        class shift_conditions_t {
            protected:
                const edge_container_t* unknowns;
                _triangulator_t& triangulator;
                
                std::vector<double> results;
                std::vector<double> partial_results;
                
                virtual void init() = 0;
            public:
                shift_conditions_t(_triangulator_t& triangulator, unsigned num_conditions) : triangulator(triangulator), results(num_conditions,0), partial_results(num_conditions,0) {}
                shift_conditions_t(_triangulator_t& triangulator) : shift_conditions_t(triangulator,3) {}
                void setUnknowns(const edge_container_t& unknowns) {
                    this->unknowns = &unknowns;
                    this->init();
                }
                virtual unsigned getNumConditions() {
                    return results.size();
                }
                virtual const edge_container_t& getConstraintEdges() = 0;
                virtual const std::vector<double>& operator()() = 0;
                virtual const std::vector<double>& partial_derivative(const edge_t& edge) {
                    // numerical differentiation, central difference quotient as default
                    // should be overloaded with analytical derivative!
                    
                    double central_len_sq = this->triangulator.getEdgeProperty(edge).length_sq;
                    double h = std::sqrt(std::numeric_limits<double>::epsilon())*central_len_sq;
                    this->triangulator.setEdgeLengthSq(edge,central_len_sq+h);
                    partial_results = this->operator ()();
                    this->triangulator.setEdgeLengthSq(edge,central_len_sq-h);
                    std::vector<double> temp = this->operator ()();
                    for(int i = 0; i<partial_results.size(); ++i) {
                        partial_results.at(i) = (partial_results.at(i) - temp.at(i))/(2*h);
                    }
                    this->triangulator.setEdgeLengthSq(edge,central_len_sq);
                    return partial_results;
                }
                
        };
        sorkin_stepper(_triangulator_t& triangulator, _regge_generator_t& regge_generator, unsigned threads, double courant) : triangulator(triangulator), regge_generator(regge_generator), threads(threads), courant(courant) { 
            if(courant > 1.0d/2.0d)
                throw std::logic_error("Courant ratio is too large!");
        }
        sorkin_stepper(_triangulator_t& triangulator, _regge_generator_t& regge_generator, unsigned threads) : sorkin_stepper(triangulator,regge_generator,threads,3) { }
        template<typename _shift_conditions_t> double step(_shift_conditions_t& shift_conditions, double dt_sq) {            
            vertex_hyper_index_pred<skeleton_t> hyper_index_pred(2,triangulator.getSkeleton());
            auto hyper_index_filtered_skeleton = triangulator.getFilteredSkeleton(hyper_index_pred);
            
            for(int i = 0; i<triangulator.getNumColors(); ++i) {
                vertex_color_pred<skeleton_t> color_pred(i,2,triangulator.getSkeleton());
                auto color_filtered_skeleton = triangulator.getFilteredSkeleton(color_pred);
                
                for(auto iter = boost::vertices(color_filtered_skeleton).first; iter != boost::vertices(color_filtered_skeleton).second; ++iter) {
                    std::vector<std::thread> thread_vec;
                    for(int j = threads; (j>0) && (iter != boost::vertices(color_filtered_skeleton).second); --j) {
                        thread_vec.push_back(std::thread(&sorkin_stepper::solve<_shift_conditions_t>,this,*iter,hyper_index_filtered_skeleton,dt_sq,shift_conditions));
                        // non threaded version:
                        //this->solve(*iter,hyper_index_filtered_skeleton,dt_sq,shift_conditions);
                        ++iter;
                    }
                    --iter;
                    for(auto thread_iter = thread_vec.begin(); thread_iter != thread_vec.end(); ++thread_iter)
                        thread_iter->join();
                }
            }
            return std::sqrt(dt_sq);
        }
        template<typename _shift_conditions_t> double step(_shift_conditions_t& shift_conditions) {
            double min_len_sq = this->getMinLength();
            double dt_sq = (courant*courant)*min_len_sq;
            return this->step(shift_conditions,dt_sq);
        }
};

#endif	/* SORKIN_STEPPER_H */


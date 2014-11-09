/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef QPL_MESH_GENERATOR_H
#define	QPL_MESH_GENERATOR_H
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include "triangulator.h"
#include "utilities.h"
#include "sorkin_stepper.h"

#include <cmath>

template<typename vertex_weight_t, typename edge_weight_t, typename triangle_weight_t, 
        typename tetrahedron_weight_t, typename pentachoron_weight_t>
class qpl_mesh_generator_t {
    private:
        typedef triangulator_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> _triangulator_t;
        typedef typename _triangulator_t::geodesics_t geodesics_t;
        typedef typename _triangulator_t::hypersurface_t hypersurface_t;
        typedef typename _triangulator_t::vertex_num_t vertex_num_t;
        typedef typename _triangulator_t::hyper_index_t hyper_index_t;
        typedef typename _triangulator_t::edge_t edge_t;
        typedef typename _triangulator_t::edge_container_t edge_container_t;
        typedef typename _triangulator_t::vertex_t vertex_t;

        
        typedef std::vector<unsigned> index_vector_t;
        typedef typename std::unordered_map<vertex_num_t,index_vector_t> index_map_t;
        typedef typename std::unordered_map<index_vector_t,vertex_num_t,vector_hash<unsigned>> reverse_index_map_t;
        
        typedef std::function<double(const unsigned,const index_vector_t&,const unsigned,const index_vector_t&)> length_sq_t;
        
        class qpl_geodesics_t : public geodesics_t {
            protected:
                const index_map_t& index_map;
                const unsigned resolution;
                length_sq_t l_sq;
            public:
                qpl_geodesics_t(const index_map_t& index_map, const unsigned resolution) : index_map(index_map), resolution(resolution) {}
                virtual double operator()(const hyper_index_t h0, const vertex_num_t index0, const hyper_index_t h1, const vertex_num_t index1) const {
                    if(std::abs(h0-h1) > 1)
                        throw std::logic_error("Trying to set edge length between Sigma0 and Sigma2.");
                    index_vector_t indices0 = index_map.at(index0);
                    index_vector_t indices1 = index_map.at(index1);
                    for(int i = 0; i<=3; ++i) {
                        int temp = std::abs(indices0.at(i) - indices1.at(i));
                        if((temp > 1) && (temp != resolution))
                            throw std::logic_error("Trying to set QPL edge with >1 coordinate difference.");
                    }
                    return l_sq(h0,indices0,h1,indices1);
                }
                void setLengthFunction(length_sq_t l_sq) {
                    this->l_sq = l_sq;
                }
        };
        
        hypersurface_t hypersurface;
        index_map_t index_map;
        reverse_index_map_t reverse_index_map;
        qpl_geodesics_t qpl_geodesics;
        
    public:
        class qpl_zeroshift_t : public sorkin_stepper<vertex_weight_t,edge_weight_t,triangle_weight_t,tetrahedron_weight_t,pentachoron_weight_t>::shift_conditions_t {
            typedef qpl_mesh_generator_t<vertex_weight_t,edge_weight_t,triangle_weight_t,tetrahedron_weight_t,pentachoron_weight_t> _qpl_mesh_generator_t;
            typedef typename sorkin_stepper<vertex_weight_t,edge_weight_t,triangle_weight_t,tetrahedron_weight_t,pentachoron_weight_t>::shift_conditions_t shift_conditions_t;
            
            protected:            
                const _qpl_mesh_generator_t& qpl;
                std::vector<std::pair<edge_t,edge_t>> zero_shift_pairs;
                edge_container_t contraint_edges;
                virtual void init() {
                    std::vector<std::pair<edge_t,unsigned>> shift_candidates;
                    for(auto iter = this->unknowns->begin(); (iter != this->unknowns->end()) && (shift_candidates.size() < 6) ; ++iter) {
                        vertex_t src = this->triangulator.getSource(*iter);
                        vertex_t tgt = this->triangulator.getTarget(*iter);
                        const index_vector_t& src_coords = qpl.getCoordinates(this->triangulator.getInputIndices(src).second);
                        const index_vector_t& tgt_coords = qpl.getCoordinates(this->triangulator.getInputIndices(tgt).second);
                        // take only 3 + 3 cubic edges
                        if(src_coords.at(0) == tgt_coords.at(0)) {
                            // compute coordiante direction
                            unsigned k = 0;
                            for(int i = 1; i<4; ++i) {
                                if(src_coords.at(i) != tgt_coords.at(i))
                                    k = i;
                            }
                            if(k == 0)
                                throw std::logic_error("QPL zeroshift condition initialized to false coordinate direction.");

                            
                            shift_candidates.push_back(std::make_pair(*iter,k));
                        }
                    }
                    if(shift_candidates.size() != 6)
                        throw std::logic_error("QPL zeroshift condition missing candidates.");
                    std::sort(shift_candidates.begin(),shift_candidates.end(),pairs_second_comparision<edge_t,unsigned>());

                    // collect them in pairs
                    for(auto iter = shift_candidates.begin(); iter != shift_candidates.end(); ++iter) {
                        edge_t e0 = iter->first;
                        ++iter;
                        edge_t e1 = iter->first;
                        zero_shift_pairs.push_back(std::make_pair(e0,e1));
                        contraint_edges.push_back(e0);
                    }
                }
            public:
                qpl_zeroshift_t(_triangulator_t& triangulator, const _qpl_mesh_generator_t& qpl) : shift_conditions_t(triangulator), qpl(qpl) {}
                qpl_zeroshift_t(const qpl_zeroshift_t& qpl_zeroshift) : qpl_zeroshift_t(qpl_zeroshift.triangulator,qpl_zeroshift.qpl) {}

                virtual const std::vector<double>& operator()() {
                    for(int i = 0; i < 3; ++i) {
                        this->results.at(i) = this->triangulator.getEdgeProperty(zero_shift_pairs.at(i).first).length_sq - this->triangulator.getEdgeProperty(zero_shift_pairs.at(i).second).length_sq;
                    }
                    return this->results;
                }
                virtual const std::vector<double>& partial_derivative(const edge_t& edge) {
                    for(int i = 0; i<3; ++i) {
                        if(zero_shift_pairs.at(i).first == edge) {
                            this->partial_results.at(i) = 1.0d;
                        }
                        else {
                            if(zero_shift_pairs.at(i).second == edge) {
                                this->partial_results.at(i) = -1.0d;
                            }
                            else {
                                this->partial_results.at(i) = 0.0d;
                            }
                        }
                        
                    }
                    return this->partial_results;
                }
                virtual const edge_container_t& getConstraintEdges() {
                    return this->contraint_edges;
                }

        };
        
        class qpl_nocond_t : public sorkin_stepper<vertex_weight_t,edge_weight_t,triangle_weight_t,tetrahedron_weight_t,pentachoron_weight_t>::shift_conditions_t {
            typedef qpl_mesh_generator_t<vertex_weight_t,edge_weight_t,triangle_weight_t,tetrahedron_weight_t,pentachoron_weight_t> _qpl_mesh_generator_t;
            typedef typename sorkin_stepper<vertex_weight_t,edge_weight_t,triangle_weight_t,tetrahedron_weight_t,pentachoron_weight_t>::shift_conditions_t shift_conditions_t;
            
            protected:            
                const _qpl_mesh_generator_t& qpl;
                edge_container_t contraint_edges;

                virtual void init() {
                    this->results.clear();
                    for(int i = 0; i<3; ++i) {
                        this->results.push_back(0);
                    }
                    this->partial_results.clear();
                    for(int i = 0; i<3; ++i) {
                        this->partial_results.push_back(0);
                    }
                }
            public:
                qpl_nocond_t(_triangulator_t& triangulator, const _qpl_mesh_generator_t& qpl) : shift_conditions_t(triangulator), qpl(qpl) {}
                qpl_nocond_t(const qpl_zeroshift_t& qpl_zeroshift) : qpl_zeroshift_t(qpl_zeroshift.triangulator,qpl_zeroshift.qpl) {}

                virtual const std::vector<double>& operator()() {
                    return this->results;
                }
                virtual const std::vector<double>& partial_derivative(const edge_t& edge) {
                    return this->partial_results;
                }
                virtual const edge_container_t& getConstraintEdges() {
                    return this->contraint_edges;
                }


        };
        
        class qpl_homogeneous_metrics_t {
                //std::function<double(const unsigned,const unsigned,const index_vector_t&)> homogeneous_length_sq;
            public:
                class metric_t {
                    protected:
                        std::vector<double> deltas;
                    
                        virtual double brace_cubic(const index_vector_t& delta_indices) const = 0;
                        virtual double brace_diagonal(const index_vector_t& delta_indices) const = 0;
                        virtual double cubic(const unsigned h, const index_vector_t& delta_indices) const = 0;
                        virtual double diagonal(const unsigned h, const index_vector_t& delta_indices) const = 0;
                    public:
                        metric_t(const double courant, const double dx, const double dy, const double dz) 
                                : deltas{courant*std::min(dx,std::min(dy,dz)),dx,dy,dz} {
                            if(courant > 1.0d/2.0d)
                                throw std::logic_error("Courant condition to larger!");
                        }
                        metric_t(const double courant, const double dx) : metric_t(courant,dx,dx,dx) { }
                        virtual double operator()(const unsigned h0, const unsigned h1, const index_vector_t& delta_indices) const {
                            bool isDiagonal = (delta_indices.at(0) != 0);
                            if(h0 != h1) {
                                bool b = true;
                                for(auto it = delta_indices.begin(); it != delta_indices.end(); ++it) {
                                    b = b && (*it == 0);
                                }
                                if(b) {
                                    return -deltas.at(0)*deltas.at(0); //timelike edge
                                }
                                else { // spacelike braces
                                    if(isDiagonal)
                                        return this->brace_diagonal(delta_indices);
                                    return this->brace_cubic(delta_indices);
                                }
                            }
                            if(isDiagonal)
                                return this->diagonal(h0,delta_indices);
                            return this->cubic(h0,delta_indices);
                        }
                        double getTimestep() {
                            return deltas.at(0);
                        }
                };
            private:
                const metric_t* metric;
            
            public:
                class kasner_metric_t : public metric_t {
                    protected:
                        const std::vector<double> p_vec;
                        
                        virtual double brace_cubic(const index_vector_t& delta_indices) const {
                            double spatial = 0;
                            double test;
                            for(int i = 1; i<4; ++i) {
                                if(delta_indices.at(i) != 0)
                                    spatial = spatial + this->deltas.at(i)*this->deltas.at(i)*(1.0d+p_vec.at(i-1)*this->deltas.at(0));
                                 
                            }
                            return spatial-this->deltas.at(0)*this->deltas.at(0);
                        }
                        virtual double brace_diagonal(const index_vector_t& delta_indices) const {
                            double spatial = 0;
                            for(int i = 1; i<4; ++i) {
                                spatial = spatial + 1.0d/4.0d*this->deltas.at(i)*this->deltas.at(i)*(1.0d+p_vec.at(i-1)*this->deltas.at(0));
                            }
                            return spatial-this->deltas.at(0)*this->deltas.at(0);
                        }
                        virtual double cubic(const unsigned h, const index_vector_t& delta_indices) const {
                            int index = 0;
                            for(int i = 1; (index == 0) && (i<4); ++i) {
                                if(delta_indices.at(i) != 0)
                                    index = i;
                            }
                            if(h == 0) {      
                                return this->deltas.at(index)*this->deltas.at(index);
                            }
                            return this->deltas.at(index)*this->deltas.at(index)*std::pow(1.0d+this->deltas.at(0),2.0d*p_vec.at(index-1));
                        }
                        virtual double diagonal(const unsigned h, const index_vector_t& delta_indices) const {
                            double temp = 0;
                            double time_dep = 0;
                            for(int i = 1; i<4; ++i) {
                                double d_sq = (this->deltas.at(i))*(this->deltas.at(i));
                                temp = temp + d_sq;
                                time_dep = time_dep + std::pow(1.0d+this->deltas.at(0),2.0d*p_vec.at(i-1))*d_sq;
                            }
                            temp = temp*1.0d/4.0d;
                            time_dep = time_dep*1.0d/4.0d;
                            if(h == 0) {
                                return temp;
                            }
                            return time_dep;//temp+time_dep;
                        }
                    public:
                        kasner_metric_t(const double courant, const double dx, const double dy, const double dz, const double px, const double py, const double pz)
                                    : metric_t(courant,dx,dy,dz), p_vec{px,py,pz} { 
                                        double min = std::min(dx,std::min(dy,dz));
                                        index_vector_t index {0,0,0,0};
                                        if(dx == min)
                                            index.at(1) = 1;
                                        else {
                                            if(dy == min)
                                                index.at(2) = 1;
                                            else {
                                                if(dz == min)
                                                    index.at(3) = 1;
                                            }
                                        }
                                        
                                        
                                        this->deltas.at(0) = courant*std::sqrt(std::min(this->diagonal(0,index),this->cubic(0,index)));   
                        }
                        kasner_metric_t(const double courant, const double dx, const double px, const double py, const double pz) 
                                    : kasner_metric_t(courant,dx,dx,dx,px,py,pz) { }
                };
                
                
                class minkowski_metric_t : public metric_t {
                    protected:
                        
                        virtual double brace_cubic(const index_vector_t& delta_indices) const {
                            double spatial = 0;
                            for(int i = 1; i<4; ++i) {
                                if(delta_indices.at(i) != 0)
                                    spatial = spatial + this->deltas.at(i)*this->deltas.at(i);
                            }
                            return -this->deltas.at(0)*this->deltas.at(0)+spatial;
                        }
                        virtual double brace_diagonal(const index_vector_t& delta_indices) const {
                            double spatial = 0;
                            for(int i = 1; i<4; ++i) {
                                spatial = spatial + this->deltas.at(i)*this->deltas.at(i);;
                            }
                            spatial = 1.0d/4.0d*spatial;
                            return -this->deltas.at(0)*this->deltas.at(0)+spatial;
                        }
                        virtual double cubic(const unsigned h, const index_vector_t& delta_indices) const {
                            int index = 0;
                            for(int i = 1; (index == 0) && (i<4); ++i) {
                                if(delta_indices.at(i) != 0)
                                    index = i;
                            }
                            return this->deltas.at(index)*this->deltas.at(index);
                        }
                        virtual double diagonal(const unsigned h, const index_vector_t& delta_indices) const {
                            double temp = 0;
                            for(int i = 1; i<4; ++i) {
                                double d_sq = (this->deltas.at(i))*(this->deltas.at(i));
                                temp = temp + d_sq;
                            }
                            temp = temp*1.0d/4.0d;
                            return temp;
                        }
                    public:
                        minkowski_metric_t(const double courant, const double dx, const double dy, const double dz)
                                    : metric_t(courant,dx,dy,dz) { }
                        minkowski_metric_t(const double courant, const double dx) 
                                    : minkowski_metric_t(courant,dx,dx,dx) { }
                };
                class LAMDA_metric_t : public metric_t {
                    protected:
                        const double a_init_sq;
                        const double H_a;
                        
                        virtual double brace_cubic(const index_vector_t& delta_indices) const {
                            double spatial = 0;
                            for(int i = 1; i<4; ++i) {
                                if(delta_indices.at(i) != 0) {
                                    spatial = spatial + a_init_sq*this->deltas.at(i)*this->deltas.at(i)*(1.0d+H_a*this->deltas.at(0));
                                }
                            }
                            return -this->deltas.at(0)*this->deltas.at(0)+spatial;
                        }
                        virtual double brace_diagonal(const index_vector_t& delta_indices) const {
                            double spatial = 0;
                            for(int i = 1; i<4; ++i) {
                                spatial = spatial + a_init_sq/4.0d*this->deltas.at(i)*this->deltas.at(i)*(1.0d+H_a*this->deltas.at(0));
                            }
                            return -this->deltas.at(0)*this->deltas.at(0)+spatial;
                        }
                        virtual double cubic(const unsigned h, const index_vector_t& delta_indices) const {
                            int index = 0;
                            for(int i = 1; (index == 0) && (i<4); ++i) {
                                if(delta_indices.at(i) != 0)
                                    index = i;
                            }
                            if(h == 0) {       
                                return a_init_sq*this->deltas.at(index)*this->deltas.at(index);
                            }
                            return a_init_sq*this->deltas.at(index)*this->deltas.at(index)*(1.0d+2.0d*H_a*this->deltas.at(0));
                        }
                        virtual double diagonal(const unsigned h, const index_vector_t& delta_indices) const {
                            double temp = 0;
                            for(auto iter = ++(this->deltas.begin()); iter != this->deltas.end(); ++iter) {
                                temp = temp + (*iter)*(*iter);
                            }
                            temp = temp*1.0d/4.0d;
                            if(h == 0) {
                                return temp;
                                return a_init_sq*temp;
                            }
                            return a_init_sq*(1.0d+2.0d*H_a*this->deltas.at(0))*temp;
                        }
                        double H_a_init(const double a_init, const unsigned resolution, const double cosmo_const) {
                            return std::sqrt(cosmo_const/3.0d);
                        }
                    public:
                        LAMDA_metric_t(const double courant, const double dx, const double dy, const double dz, 
                                const double a_init, const double cosmo_const, const unsigned resolution) 
                                    : metric_t(courant,dx,dy,dz), a_init_sq(a_init*a_init), H_a(this->H_a_init(a_init,resolution,cosmo_const)) { 
                                        double min = std::min(dx,std::min(dy,dz));
                                        index_vector_t index {0,0,0,0};
                                        if(dx == min)
                                            index.at(1) = 1;
                                        else {
                                            if(dy == min)
                                                index.at(2) = 1;
                                            else {
                                                if(dz == min)
                                                    index.at(3) = 1;
                                            }
                                        }
                                        
                                        this->deltas.at(0) = courant*std::sqrt(std::min(this->diagonal(0,index),this->cubic(0,index)));
                                    }
                        LAMDA_metric_t(const double courant, const double dx, 
                                const double a_init, const double cosmo_const, const unsigned resolution)
                                    : LAMDA_metric_t(courant,dx,dx,dx,a_init,cosmo_const,resolution) { }

                };
                
                qpl_homogeneous_metrics_t(const metric_t& metric) : metric(&metric) { }
                virtual double operator()(const hyper_index_t h0, const index_vector_t indices0, const hyper_index_t h1, const index_vector_t indices1) const {
                    index_vector_t delta_indices;
                    for(int i = 0; i<=3; ++i) {
                        unsigned temp = std::abs(indices0.at(i) - indices1.at(i));
                        if(temp > 1)
                            temp = 1;
                        delta_indices.push_back(temp);
                    }
                    return (*metric)(h0,h1,delta_indices);
                }
        };
        
        qpl_mesh_generator_t(unsigned resolution) : qpl_geodesics(index_map, resolution) {
            typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_t;
            if(resolution < 3)
                throw std::logic_error("QPL resolution to low, leading to falsely identified simplices.");
            graph_t skeleton;
            for(unsigned i = 0; i <= resolution; ++i) {
                for(unsigned j = 0; j <= resolution; ++j) {
                    for(unsigned k = 0; k <= resolution; ++k) {
                        vertex_num_t v0 = boost::add_vertex(skeleton);
                        vertex_num_t v1 = boost::add_vertex(skeleton);
                        
                        std::vector<unsigned> vec0;
                        vec0.push_back(0);
                        vec0.push_back(i);
                        vec0.push_back(j);
                        vec0.push_back(k);
                        std::vector<unsigned> vec1;
                        vec1.push_back(1);
                        vec1.push_back(i);
                        vec1.push_back(j);
                        vec1.push_back(k);
                        
                        if(!index_map.emplace(std::make_pair(v0,vec0)).second)
                            throw std::logic_error("Something went wrong while inserting the internal QPL index.");
                        if(!index_map.emplace(std::make_pair(v1,vec1)).second)
                            throw std::logic_error("Something went wrong while inserting the internal QPL index.");
                        if(!reverse_index_map.emplace(std::make_pair(vec0,v0)).second)
                            throw std::logic_error("Something went wrong while inserting the internal QPL index.");
                        if(!reverse_index_map.emplace(std::make_pair(vec1,v1)).second)
                            throw std::logic_error("Something went wrong while inserting the internal QPL index.");

                    }
                }
            }
            for(auto iter = reverse_index_map.begin(); iter != reverse_index_map.end(); ++iter) {
                std::vector<unsigned> temp_coords = iter->first;
                for(auto it = ++temp_coords.begin(); it != temp_coords.end(); ++it) {
                    unsigned temp = *it;
                    if(temp == resolution)
                        *it = 0;
                    else 
                        *it = temp + 1;
                    vertex_num_t index = reverse_index_map.at(temp_coords);
                    if(boost::edge(iter->second,index,skeleton).second)
                        throw std::logic_error("Trying to add an already existent edge to the QPL.");
                    if(!boost::add_edge(iter->second,index,skeleton).second)
                        throw std::logic_error("Adding an edge to the QPL failed.");
                    *it = temp;
                }
                if(iter->first.front() == 1) {
                    for(unsigned i = 0; i <= 1; ++i) {
                        for(unsigned j = 0; j <= 1; ++j) {
                            for(unsigned k = 0; k <= 1; ++k) {
                                std::vector<unsigned> diag_coords = temp_coords;
                                diag_coords.at(0) = 0;
                                diag_coords.at(1) = diag_coords.at(1)+i;
                                diag_coords.at(2) = diag_coords.at(2)+j;
                                diag_coords.at(3) = diag_coords.at(3)+k;
                                for(int l = 1; l < 4; ++l) {
                                    if(diag_coords.at(l) > resolution)
                                        diag_coords.at(l) = 0;
                                }
                                vertex_num_t index = reverse_index_map.at(diag_coords);
                                if(boost::edge(iter->second,index,skeleton).second)
                                    throw std::logic_error("Trying to add an already existent edge to the QPL.");
                                if(!boost::add_edge(iter->second,index,skeleton).second)
                                    throw std::logic_error("Adding an edge to the QPL failed.");
                            }
                        }
                    }
                }
            }
            //std::cout << "Num Verts " << boost::num_vertices(skeleton) << " Num Edges " << boost::num_edges(skeleton) 
            //        << " E/V = " << boost::num_edges(skeleton)/boost::num_vertices(skeleton) << std::endl;
            clique_collector_visitor<graph_t> clique_collector;
            boost::bron_kerbosch_all_cliques(skeleton,clique_collector.getVisitor());
            for(auto iter = clique_collector.getCliques().begin(); iter != clique_collector.getCliques().end(); ++iter) {
                hypersurface.addTetrahedron(*iter);
            }
        }
        hypersurface_t& getHypersurface() {
            return hypersurface;
        }
        
        const geodesics_t& getGeodesics(length_sq_t length_sq) {
            qpl_geodesics.setLengthFunction(length_sq);
            return qpl_geodesics;
        }
        const index_vector_t& getCoordinates(vertex_num_t index) const {
            typename index_map_t::const_iterator it = index_map.find(index);
            if(it == index_map.end())
                throw std::logic_error("Could not find QPL coordinates.");
            return it->second;
        }
};

#endif	/* QPL_MESH_GENERATOR_H */


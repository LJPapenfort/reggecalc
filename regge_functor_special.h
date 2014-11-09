/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef REGGE_FUNCTOR_SPECIAL_H
#define	REGGE_FUNCTOR_SPECIAL_H

#include <vector>
#include <limits>
#include "regge_generator.h"
#include "utilities.h"

template<typename vertex_weight_t = empty_weight, typename edge_weight_t = empty_weight, typename triangle_weight_t = empty_weight, 
        typename tetrahedron_weight_t = empty_weight, typename pentachoron_weight_t = empty_weight>
class regge_term_t : public regge_functor_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> {
    private:
        typedef regge_functor_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> _regge_functor_t;
        typedef typename _regge_functor_t::edge_property_t edge_property_t;
        typedef typename _regge_functor_t::edge_t edge_t;
        typedef typename _regge_functor_t::edge_container_t edge_container_t;
        typedef typename _regge_functor_t::_triangulator_t::triangle_t triangle_t;
        typedef typename _regge_functor_t::_triangulator_t::pentachoron_t pentachoron_t;
        typedef typename _regge_functor_t::Matrix4d Matrix4d;
    protected:
        inline double getm1m1(const Matrix4d& metric) const {
            int temp = signum(metric(2,2)-metric(2,3)*metric(2,3)/metric(3,3));
            return temp;
        }
        inline double getm1m2(const Matrix4d& metric) const {
            int temp = signum(metric(2,2)*metric(3,3)-metric(2,3)*metric(2,3));
            return -temp*metric(2,3)/std::sqrt(std::abs(metric(2,2)*metric(3,3)));
        }
        inline double getn1m2(const Matrix4d& metric) const {
            int temp = signum(metric(3,3));
            return -temp*std::sqrt(std::abs(1-metric(2,3)*metric(2,3)/(metric(2,2)*metric(3,3))));
        }
        double getTimelikeDAngle(const Matrix4d& metric) const {
            return std::acos(this->getm1m2(metric));
        }
        double getSpacelikeDAngle(const Matrix4d& metric) const {
            double m1m2 = this->getm1m2(metric);
            double n1m2 = this->getn1m2(metric);
            double m1m1 = this->getm1m1(metric);
            double rho;
            if(std::abs(m1m2) < std::abs(n1m2)) {
                rho = signum(n1m2)*m1m2;
            }
            else {
                rho = signum(m1m2)*n1m2;
            }
            return m1m1*std::asinh(rho);
        }
        double getTriangleDeterminant(const triangle_t& triangle) const {
            double temp = 0.0d;
            for(auto iter = triangle.edges.begin(); iter != triangle.edges.end(); ++iter) {
                auto next = std::next(iter);
                if(next == triangle.edges.end())
                    next = triangle.edges.begin();
                double l1_sq = this->triangulator.getEdgeProperty(*iter).length_sq;
                double l2_sq = this->triangulator.getEdgeProperty(*next).length_sq;
                temp += 2.0d*l1_sq*l2_sq - l1_sq*l1_sq;
            }
            return 1.0d/4.0d*temp;
        }

    public:
        virtual double operator()(const edge_t& edge) const {
            const edge_property_t& e_prop = this->triangulator.getEdgeProperty(edge);
            std::vector<double> results;
            // calculate regge term for every hinge (edge,x,y)
            for(auto iter = e_prop.triangles.begin(); iter != e_prop.triangles.end(); ++iter) {
                const triangle_t& triangle = iter->get();
                // check if triangle is space/timelike and collect length for area term
                bool isTimelike = false;
                double temp_area_term = 0.0d;
                for(auto tri_iter = triangle.edges.begin(); tri_iter != triangle.edges.end(); ++tri_iter) {
                    double length_sq = this->triangulator.getEdgeProperty(*tri_iter).length_sq;
                    if(length_sq < 0)
                        isTimelike = true;
                    if(*tri_iter == edge)
                        temp_area_term += -length_sq;
                    else
                        temp_area_term += length_sq;
                }
                // calculate defect
                double defect;
                if(isTimelike) {
                    defect = 2.0d*this->PI();
                }
                else {
                    defect = 0.0d;
                }

                // get dihedral angle from standard frame of pentachora
                for(auto in_iter = triangle.pentachorons.begin(); in_iter != triangle.pentachorons.end(); ++in_iter) {
                    const pentachoron_t& penta = in_iter->get();
                    // get and transform metric with triangle to origin
                    
                    Matrix4d metric = this->getMetric(penta.edges);
                    
                    if(metric.hasNaN())
                        throw std::logic_error("Metric got NaN");
                    
                    this->transformMetric(metric,penta.edges,triangle);
                    if(metric.hasNaN())
                        throw std::logic_error("Metric transformed got NaN");
                    
                    // numerical inversion of 4x4 matrix fine with Eigen
                    metric = metric.inverse().eval();
                    if(metric.hasNaN())
                        throw std::logic_error("Metric inverse got NaN");
                    
                    
                    double phi;
                    if(isTimelike) {
                        phi = this->getTimelikeDAngle(metric);
                    }
                    else {
                        phi = this->getSpacelikeDAngle(metric);

                    }
                    defect = defect - phi;
                }
                // calculate triangle area term
                double tri_det = this->getTriangleDeterminant(triangle);
                double tri_area = 1.0d/2.0d*std::sqrt(std::abs(tri_det));
                double area_term = 1.0d/(16.0d*tri_area)*signum(tri_det)*temp_area_term;
                
                // put result in vector
                results.push_back(defect*area_term);
            }
            return this->sumUp(results);
        }
        using _regge_functor_t::regge_functor_t;
};

template<typename vertex_weight_t = empty_weight, typename edge_weight_t = empty_weight, typename triangle_weight_t = empty_weight, 
        typename tetrahedron_weight_t = empty_weight, typename pentachoron_weight_t = empty_weight>
class cosmological_const_term_t : public regge_functor_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> {
    private:
        typedef regge_functor_t<vertex_weight_t, edge_weight_t, triangle_weight_t, tetrahedron_weight_t, pentachoron_weight_t> _regge_functor_t;
        typedef typename _regge_functor_t::edge_property_t edge_property_t;
        typedef typename _regge_functor_t::edge_t edge_t;
        typedef typename _regge_functor_t::edge_container_t edge_container_t;
        typedef typename _regge_functor_t::_triangulator_t _triangulator_t;
        typedef typename _triangulator_t::pentachoron_t pentachoron_t;
        typedef typename _regge_functor_t::Matrix4d Matrix4d;
        
        double cosmo_const;
        
        double computeVolume(const double determinant) const {
            return 1.0d/24.0d*std::sqrt(std::abs(determinant));
        }
        
        double computeTerm(const edge_t& edge, const pentachoron_t& p) const {
            Matrix4d metric = this->getMetric(p.edges);
            Matrix4d metric_deriv;
            
            
            auto index_pair = this->getEdgeIndices(edge,p.edges);
            unsigned k = index_pair.first;
            unsigned l = index_pair.second;
            
            for(int i = 0; i < 4; ++i) {
                for(int j = 0; j < 4; ++j) {
                    metric_deriv(i,j) = this->getdMetricComponent(i+1,j+1,k,l);
                }
            }
            
            Matrix4d temp = metric.inverse()*metric_deriv; 
            return this->computeVolume(metric.determinant())*temp.trace();  
        }
        double getdMetricComponent(unsigned i, unsigned j, unsigned k, unsigned l) const {
            // create mapping e -> ij
            double temp = 0;
            if((i != 0) && (l == i) && (k == 0))
                temp = temp + 1.0d;
            if((j != 0) && (l == j) && (k == 0))
                temp = temp + 1.0d;
            if(((i == k) && (j == l)) || ((i == l) && (j == k)))
                temp = temp - 1.0d;
            return 1.0d/2.0d*temp;
        }
        std::pair<unsigned,unsigned> getEdgeIndices(const edge_t& edge, const edge_container_t& edges) const { 
            unsigned i = 0;
            unsigned j = 1;
            for(auto iter = edges.begin(); iter != edges.end(); ++iter) {
                if(*iter == edge)
                    iter = --edges.end();
                else {
                    if(j < 4)
                        j++;
                    else {
                        i++;
                        j = i+1;
                    }
                }
            }
            return std::make_pair(i,j);
        }
    public:
        cosmological_const_term_t(_triangulator_t& t, double cosmological_constant) : _regge_functor_t(t), cosmo_const(cosmological_constant) {}
        cosmological_const_term_t(const cosmological_const_term_t& cosmological_const_term) : _regge_functor_t(cosmological_const_term.triangulator), cosmo_const(cosmological_const_term.cosmo_const) {}
        virtual double operator()(const edge_t& edge) const {
            const edge_property_t& edge_property = this->triangulator.getEdgeProperty(edge);
            std::vector<double> results;
            
            for(auto iter = edge_property.pentachorons.begin(); iter != edge_property.pentachorons.end(); ++iter) {
                double term = this->computeTerm(edge,iter->get());
                results.push_back(term);
            }
            return -cosmo_const/2.0d*this->sumUp(results);
        }
};

#endif	/* REGGE_FUNCTOR_SPECIAL_H */


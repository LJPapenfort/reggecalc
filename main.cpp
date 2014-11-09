/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */


#include "triangulator.h"
#include <iostream>
#include "utilities.h"
#include "qpl_mesh_generator.h"
#include "sorkin_stepper.h"


int main(int argc, char** argv) {
    
    typedef triangulator_t<empty_weight,empty_weight,empty_weight,empty_weight,empty_weight> _triangulator_t;
    typedef typename _triangulator_t::edge_t edge_t;
    typedef typename _triangulator_t::edge_container_t edge_container_t;
    typedef typename _triangulator_t::vertex_t vertex_t;
    typedef regge_generator_t<empty_weight,empty_weight,empty_weight,empty_weight,empty_weight> _regge_generator_t;
    typedef regge_term_t<empty_weight,empty_weight,empty_weight,empty_weight,empty_weight> _regge_term_t;
    typedef cosmological_const_term_t<empty_weight,empty_weight,empty_weight,empty_weight,empty_weight> _cosmo_term_t;
    
    double cosmo_const = 1.0d;  
    double a_init = 1.0d;
    unsigned resolution = 3;
    
    double courant = 1.0d/2.5d;
    double dx = 0.05d;
    
    typedef qpl_mesh_generator_t<empty_weight,empty_weight,empty_weight,empty_weight,empty_weight> qpl_t;
    qpl_t qpl(resolution);
    
    qpl_t::qpl_homogeneous_metrics_t::kasner_metric_t metric(courant,dx,2.0d/3.0d,2.0d/3.0d,-1.0d/3.0d);
    //qpl_t::qpl_homogeneous_metrics_t::kasner_metric_t metric(courant,dx,0,0,1);
    //qpl_t::qpl_homogeneous_metrics_t::kasner_metric_t metric(courant,dx,1.0d/2.0d,1.0d/4.0d*(1.0d-std::sqrt(5.0d)),1.0d/4.0d*(1.0d+std::sqrt(5.0d)));
    //qpl_t::qpl_homogeneous_metrics_t::kasner_metric_t metric(courant,dx,1.0d/4.0d,1.0d/8.0d*(3.0d+std::sqrt(21.0d)),1.0d/8.0d*(3.0d-std::sqrt(21.0d)));
    //qpl_t::qpl_homogeneous_metrics_t::kasner_metric_t metric(courant,dx,-1.0d/4.0d,1.0d/8.0d*(5.0d-std::sqrt(5.0d)),1.0d/8.0d*(5.0d+std::sqrt(5.0d)));
    //qpl_t::qpl_homogeneous_metrics_t::kasner_metric_t metric(courant,dx,-1.0d/5.0d,1.0d/5.0d*(3.0d-std::sqrt(3.0d)),1.0d/5.0d*(3.0d+std::sqrt(3.0d)));
    //qpl_t::qpl_homogeneous_metrics_t::kasner_metric_t metric(courant,dx,1.0d/5.0d,2.0d/5.0d*(1.0d-std::sqrt(2.0d)),2.0d/5.0d*(1.0d+std::sqrt(2.0d)));

    //qpl_t::qpl_homogeneous_metrics_t::LAMDA_metric_t metric(courant,dx,a_init,cosmo_const,resolution);
    //qpl_t::qpl_homogeneous_metrics_t::minkowski_metric_t metric(courant,dx);

    qpl_t::qpl_homogeneous_metrics_t qpl_metric(metric);
    _triangulator_t TR(qpl.getHypersurface(),qpl.getGeodesics(qpl_metric));
    
    qpl_t::qpl_zeroshift_t zero_shift(TR,qpl);
    
    _regge_term_t r_term(TR);
    _cosmo_term_t c_term(TR,cosmo_const); 
    _regge_generator_t gen(TR,r_term); //,c_term);
    

    /*
    for(auto iter = TR.getEdges().first; iter != TR.getEdges().second; ++iter) {
        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        if((TR.getVertexProperty(src).hyper_index == 0) && (TR.getVertexProperty(tgt).hyper_index == 1)) {
            if((TR.getVertexProperty(src).next == tgt) && (TR.getVertexProperty(tgt).last == src)) {
                std::cout << "Result for " << *iter << ": " << gen(*iter) << std::endl;
                std::cout << "Length sq for " << *iter << ": " << TR.getEdgeProperty(*iter).length_sq << std::endl;
                std::cout << "Color: " << TR.getVertexProperty(src).color << std::endl;
                std::cout << "Num Penta: " << TR.getEdgeProperty(*iter).triangles.size() << std::endl;
            }
        }
        
    }
    */
    
    edge_container_t cubics;
    for(auto iter = TR.getEdges().first; iter != TR.getEdges().second; ++iter ){

        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        auto src_pair = TR.getInputIndices(src);
        auto tgt_pair = TR.getInputIndices(tgt);

        auto src_qpl = qpl.getCoordinates(src_pair.second);
        auto tgt_qpl = qpl.getCoordinates(tgt_pair.second);
        
        if((src_qpl.at(0) == tgt_qpl.at(0)) && (src_pair.first == tgt_pair.first) && (src_pair.first == 1)) {
            cubics.push_back(*iter);
        }
    }
    
    edge_container_t diagonals;
    for(auto iter = TR.getEdges().first; iter != TR.getEdges().second; ++iter ){

        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        auto src_pair = TR.getInputIndices(src);
        auto tgt_pair = TR.getInputIndices(tgt);

        auto src_qpl = qpl.getCoordinates(src_pair.second);
        auto tgt_qpl = qpl.getCoordinates(tgt_pair.second);
        
        if((src_qpl.at(0) != tgt_qpl.at(0)) && (src_pair.first == tgt_pair.first) && (src_pair.first == 1)) {
            diagonals.push_back(*iter);
        }
    }
    
    edge_container_t diagonal_braces;
    for(auto iter = TR.getEdges().first; iter != TR.getEdges().second; ++iter ){

        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        auto src_pair = TR.getInputIndices(src);
        auto tgt_pair = TR.getInputIndices(tgt);

        auto src_qpl = qpl.getCoordinates(src_pair.second);
        auto tgt_qpl = qpl.getCoordinates(tgt_pair.second);
        
        if((src_qpl.at(0) != tgt_qpl.at(0)) && (src_pair.first != tgt_pair.first) && (src_pair.first != 0) && (TR.getVertexProperty(tgt).last != src)) {
            diagonal_braces.push_back(*iter);
        }
    }
    
    edge_container_t cubic_braces;
    for(auto iter = TR.getEdges().first; iter != TR.getEdges().second; ++iter ){

        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        auto src_pair = TR.getInputIndices(src);
        auto tgt_pair = TR.getInputIndices(tgt);

        auto src_qpl = qpl.getCoordinates(src_pair.second);
        auto tgt_qpl = qpl.getCoordinates(tgt_pair.second);
        
        if((src_qpl.at(0) == tgt_qpl.at(0)) && (src_pair.first != tgt_pair.first) && (src_pair.first != 0) && (TR.getVertexProperty(tgt).last != src)) {
            cubic_braces.push_back(*iter);
        }
    }
    
    edge_container_t timelikes;
    for(auto iter = TR.getEdges().first; iter != TR.getEdges().second; ++iter ){

        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        auto src_pair = TR.getInputIndices(src);
        auto tgt_pair = TR.getInputIndices(tgt);

        auto src_qpl = qpl.getCoordinates(src_pair.second);
        auto tgt_qpl = qpl.getCoordinates(tgt_pair.second);
        
        if((tgt_pair.first == 2) && (TR.getVertexProperty(tgt).last == src)) {
            timelikes.push_back(*iter);
        }
    }
    
    edge_container_t x_cubics;
    edge_container_t y_cubics;
    edge_container_t z_cubics;
    for(auto iter = cubics.begin(); iter != cubics.end(); ++iter) {
        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        auto src_pair = TR.getInputIndices(src);
        auto tgt_pair = TR.getInputIndices(tgt);

        auto src_qpl = qpl.getCoordinates(src_pair.second);
        auto tgt_qpl = qpl.getCoordinates(tgt_pair.second);
        
        if(src_qpl.at(1) != tgt_qpl.at(1))
            x_cubics.push_back(*iter);
        if(src_qpl.at(2) != tgt_qpl.at(2))
            y_cubics.push_back(*iter);
        if(src_qpl.at(3) != tgt_qpl.at(3))
            z_cubics.push_back(*iter);
    }
    
    edge_container_t x_cubic_braces;
    edge_container_t y_cubic_braces;
    edge_container_t z_cubic_braces;
    for(auto iter = cubics.begin(); iter != cubics.end(); ++iter) {
        vertex_t src = TR.getSource(*iter);
        vertex_t tgt = TR.getTarget(*iter);
        auto src_pair = TR.getInputIndices(src);
        auto tgt_pair = TR.getInputIndices(tgt);

        auto src_qpl = qpl.getCoordinates(src_pair.second);
        auto tgt_qpl = qpl.getCoordinates(tgt_pair.second);
        
        if(src_qpl.at(1) != tgt_qpl.at(1))
            x_cubic_braces.push_back(*iter);
        if(src_qpl.at(2) != tgt_qpl.at(2))
            y_cubic_braces.push_back(*iter);
        if(src_qpl.at(3) != tgt_qpl.at(3))
            z_cubic_braces.push_back(*iter);
    }

    std::cout << "Got " << std::endl
            << cubics.size() << " cubics" << std::endl
            << x_cubics.size() << " x cubics" << std::endl
            << y_cubics.size() << " y cubics" << std::endl
            << z_cubics.size() << " z cubics" << std::endl
            << diagonals.size() << " diagonals" << std::endl
            << cubic_braces.size() << " cubic braces" << std::endl
            << x_cubic_braces.size() << " x cubic braces" << std::endl
            << y_cubic_braces.size() << " y cubic braces" << std::endl
            << z_cubic_braces.size() << " z cubic braces" << std::endl
            << diagonal_braces.size() << " diagonal braces" << std::endl
            << timelikes.size() << " timelikes" << std::endl
            << TR.getNumVertices() << " vertices" << std::endl
            << TR.getNumEdges() << " edges" << std::endl
            << TR.getNumTriangles() << " triangles" << std::endl
            << TR.getNumTetrahedra() << " tetrahedra" << std::endl
            << TR.getNumPentachora() << " pentachora" << std::endl;

    
    edge_container_t obs_edges;
    obs_edges.push_back(x_cubics.front());
    obs_edges.push_back(y_cubics.front());
    obs_edges.push_back(z_cubics.front());
    
    obs_edges.push_back(diagonals.front());
    
    obs_edges.push_back(x_cubic_braces.front());
    obs_edges.push_back(y_cubic_braces.front());
    obs_edges.push_back(z_cubic_braces.front());
    
    obs_edges.push_back(diagonal_braces.front());
    obs_edges.push_back(timelikes.front());

    sorkin_stepper<empty_weight,empty_weight,empty_weight,empty_weight,empty_weight> sorkin(TR,gen,8,courant);
    double dt = sorkin.step(zero_shift);    
    
    double t = 1.0d+std::sqrt(-TR.getEdgeProperty(TR.getEdge(TR.getVertex(0,0),TR.getVertex(1,0)).first).length_sq);
    std::cout << t << ",";
    for(auto iter = obs_edges.begin(); iter != obs_edges.end(); ++iter) {
        std::cout << TR.getEdgeProperty(*iter).length_sq;
        if(std::next(iter) != obs_edges.end())
            std::cout << ",";
    }
    std::cout << std::endl;
    
    TR.swapHypersurfaces();
    
    for(int j = 0; j < 10000; ++j) {        
        t = t + dt;
        dt = sorkin.step(zero_shift); //,dt*dt);
        
        std::cout << t << ",";
        for(auto iter = obs_edges.begin(); iter != obs_edges.end(); ++iter) {
            std::cout << TR.getEdgeProperty(*iter).length_sq;
            if(std::next(iter) != obs_edges.end())
                std::cout << ",";
        }
        std::cout << std::endl;
        
        TR.swapHypersurfaces();        
    }
    return 0;
}


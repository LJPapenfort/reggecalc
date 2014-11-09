/* 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */

#ifndef UTILITIES_H
#define	UTILITIES_H
#include <unordered_set>

// helpers


template<typename T> int signum(T t) {
    int temp = (T(0) < t) - (t < T(0));
    return temp;
}


unsigned kronecker_delta(int i, int j) {
    return i==j;
}

template<typename T1, typename T2> struct pairs_second_comparision {
        bool operator()(const std::pair<T1,T2> pair1, const std::pair<T1,T2> pair2) {
            return pair1.second < pair2.second;
        }
};

// hashes

template<typename T, typename hash_t = std::hash<T>>
class vector_hash {
    typedef std::vector<T> vector_t;
    typedef typename vector_t::const_iterator c_iterator;
    public:
        std::size_t operator()(const vector_t& vec) const {
            std::size_t seed = 0;
            for(c_iterator iter = vec.begin(); iter != vec.end(); ++iter) {
                boost::hash_combine(seed, hash_t()(*iter));
            }
            return seed;
        }
};

template<typename T, typename S>
class pair_hash {
    typedef std::pair<T,S> pair_t;
    public:
        std::size_t operator()(const pair_t& pair) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, pair.first);
            boost::hash_combine(seed, pair.second);
            return seed;
        }
};

template<typename skeleton_t>
class edge_t_hash {
    private:
        typedef typename boost::graph_traits<skeleton_t>::edge_descriptor edge_t;
    public:
        std::size_t operator()(const edge_t& e) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, std::min(e.m_source,e.m_target));
            boost::hash_combine(seed, std::max(e.m_source,e.m_target));
            return seed;
        }
};

// predicates

template<typename graph_t>
class simplex_swap_pred {
    typedef boost::graph_traits<graph_t> traits;
    typedef typename traits::edge_descriptor edge_t;
    typedef typename traits::vertex_descriptor vertex_t;
    
    const graph_t& graph;
    public:
        simplex_swap_pred(const graph_t& graph) : graph(graph) {}
        bool operator()(const edge_t& edge) const {
            if(graph[boost::source(edge,graph)].hyper_index < 1)
                return true;
            if(graph[boost::target(edge,graph)].hyper_index < 1)
                return true;
            return false;
        }
};

template<typename graph_t>
class adjacent_vertex_pred {
    private:
        typedef boost::graph_traits<graph_t> traits;
        typedef typename traits::vertex_descriptor v_descriptor;
        typedef typename traits::adjacency_iterator a_iterator;
        std::unordered_set<v_descriptor> v_set;
    public:
        bool operator()(const v_descriptor& v) const {
            return v_set.find(v) != v_set.end();
        }
        adjacent_vertex_pred() { }
        adjacent_vertex_pred(const v_descriptor& v, const graph_t& g) {
            a_iterator a_iter, a_end;
            for(boost::tie(a_iter,a_end) = boost::adjacent_vertices(v,g); a_iter != a_end; ++a_iter) {
                v_descriptor a_vertex = *a_iter;
                bool b = (!v_set.insert(*a_iter).second);
                    //throw something with b
                
            }
        }
        typename std::unordered_set<v_descriptor>::size_type num_adjacent_vertices() const {
            return v_set.size();
        }        
};

template<typename graph_t>
class vertex_color_pred {
    private:
        typedef boost::graph_traits<graph_t> traits;
        typedef typename traits::vertex_descriptor v_descriptor;
        typedef typename traits::vertices_size_type color_t;
        const graph_t* graph;
        color_t color;
        unsigned short index;
    public:
        vertex_color_pred() : color(0) , index(0) { }
        vertex_color_pred(color_t c, unsigned short i, const graph_t& g) : color(c), index(i), graph(&g) { }
        bool operator()(const v_descriptor& v) const {
            return (index == (*graph)[v].hyper_index) && (color == (*graph)[v].color);
        } 
        
        
};

template<typename graph_t>
class vertex_hyper_index_pred {
    private:
        typedef boost::graph_traits<graph_t> traits;
        typedef typename traits::vertex_descriptor v_descriptor;
        typedef typename traits::vertices_size_type color_t;
        const graph_t* graph;
        unsigned short index;
    public:
        vertex_hyper_index_pred() :  index(0) { }
        vertex_hyper_index_pred(unsigned short h, const graph_t& g) : index(h), graph(&g) { }
        bool operator()(const v_descriptor& v) const {
            return (index == (*graph)[v].hyper_index);
        } 
        
        
};

template<typename graph_t>
class vertex_vector_pred {
    private:
        typedef boost::graph_traits<graph_t> traits;
        typedef typename traits::vertex_descriptor v_descriptor;
        typedef typename traits::vertices_size_type color_t;
        const std::vector<v_descriptor>* vertex_vector;
    public:
        vertex_vector_pred() { }
        vertex_vector_pred(const std::vector<v_descriptor>& vert_vec) : vertex_vector(&vert_vec) { }
        bool operator()(const v_descriptor& v) const {
            return std::find((*vertex_vector).begin(),(*vertex_vector).end(),v) != (*vertex_vector).end();
        }
};

// visitors
template<typename graph_t>
class clique_collector_visitor {
        typedef typename boost::graph_traits<graph_t>::vertex_descriptor v_descriptor;
        typedef std::vector<v_descriptor> vertex_container_t;
        typedef std::vector<vertex_container_t> clique_container_t;
        clique_container_t clique_container;
        class bk_visitor {
            clique_container_t& cc;
            public:
                bk_visitor(clique_container_t& input_cc) : cc(input_cc) {}
                template<typename clique_t> void clique(const clique_t& c, const graph_t& g) {
                    vertex_container_t vertex_container;
                    for(auto iter = c.begin(); iter != c.end(); ++iter) {
                        vertex_container.push_back(*iter);
                    }
                    cc.push_back(std::move(vertex_container));
                }
        };
        bk_visitor vis;
    public:
        clique_collector_visitor() : vis(clique_container) {}
        clique_container_t& getCliques() {
            return clique_container;
        }
        bk_visitor& getVisitor() {
            return vis;
        }
};

class clique_counter_visitor {
    int cliques_num;
    int maximum_clique_num;
    int minimum_clique_num;
    public:
        clique_counter_visitor() : cliques_num(0), maximum_clique_num(0), minimum_clique_num(0) { }
        template <typename clique_t, typename graph_t>
        void clique(const clique_t& c, const graph_t& g) {
            typedef typename clique_t::const_iterator iterator;
            iterator v_iter, v_end = c.end();
            int i = 0;
            for(v_iter = c.begin(); v_iter != v_end; ++v_iter)
            {
                i++;
            }
            if(i>maximum_clique_num) {
                maximum_clique_num = i;
                if(minimum_clique_num == 0)
                    minimum_clique_num = i;
            }
            if(i<minimum_clique_num)
                minimum_clique_num = i;
            cliques_num++;
            
        }
        int cliques_number() {
            return cliques_num;
        }
        int maximum_clique_number() {
            return maximum_clique_num;
        }
        int minimum_clique_number() {
            return minimum_clique_num;
        }
};

#endif	/* UTILITIES_H */


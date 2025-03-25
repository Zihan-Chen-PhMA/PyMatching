#ifndef DIJKSTRA_GRAPH_H
#define DIJKSTRA_GRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include "pymatching/sparse_blossom/flooder/detector_node.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"

namespace dijkstra{


    struct VertexData {
        std::string v_tag;
        VertexData();
        VertexData(std::string v_tag);    
    };

    struct WeightedEdgeData {
        int64_t weight;
        int64_t weight_copy;
        int64_t *flooded_weight_st;
        int64_t *flooded_weight_ts;
        WeightedEdgeData();
        WeightedEdgeData(int64_t weight);
    };

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
            VertexData, WeightedEdgeData> flooded_graph_t;
    typedef boost::graph_traits<flooded_graph_t>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<flooded_graph_t>::edge_descriptor edge_descriptor;
    typedef boost::graph_traits<flooded_graph_t>::edge_iterator edge_iterator;
    typedef boost::graph_traits<flooded_graph_t>::vertex_iterator vertex_iterator;
    typedef boost::graph_traits<flooded_graph_t>::out_edge_iterator out_edge_iterator;
    typedef boost::property_map<flooded_graph_t, boost::vertex_index_t>::type pmap_v_index;



    class DijkstraNodeVisitor : public boost::default_dijkstra_visitor  
    {
        public:
        DijkstraNodeVisitor(const vertex_descriptor& target,
                        std::vector<vertex_descriptor>& discovered_nodes);

        void examine_vertex(vertex_descriptor v, const flooded_graph_t &g)
        {
            if (v == target){
                throw (2);
            }   
        }
        void discover_vertex(const vertex_descriptor& v, const flooded_graph_t &g){
            discovered_nodes.push_back(v);
        }
        vertex_descriptor target;
        std::vector<vertex_descriptor>& discovered_nodes;
    };

    class Node {
        public:
        int64_t local_radius;
        std::vector<int64_t> neighbor_weights;
        std::vector<int64_t> neighbor_flooded_weights;
        Node() : local_radius(0) {}
    };


    class SoftOutputDijkstra {
        public:
        flooded_graph_t fl_matching_graph;
        std::vector<pm::DetectorNode*> flooded_nodes;
        std::vector<size_t> flooded_node_indices; 
        std::vector<std::pair<size_t,size_t>> cycle_endpoint_pairs;
        std::vector<std::pair<size_t,size_t>> cycle_endpoint_pairs_mono;
        std::vector<Node> nodes;
        
        

        SoftOutputDijkstra();
        SoftOutputDijkstra(pm::Mwpm &mwpm);

        flooded_graph_t mwpm_to_dijkstra_graph(pm::Mwpm &mwpm);

        void add_edge_from_mwpm_to_graph(size_t source_index, size_t target_index, pm::Mwpm &mwpm,
                                         flooded_graph_t &graph);

        void add_boundary_node(size_t boundary_index);

        void add_boundary_edge_from_mwpm(size_t inner_index, size_t boundary_index, pm::Mwpm &mwpm);

        void add_image_node(size_t original_index, size_t image_index); 

        void redirect_edge_to_image(size_t original_index, size_t image_index, size_t target_index);

        void add_boundary_edge_to_image_from_mwpm(size_t original_index, size_t image_index, size_t boundary_index, pm::Mwpm &mwpm);

        void copy_edge_to_image(size_t original_s_index, size_t original_t_index, size_t image_s_index, size_t image_t_index);

        void add_cycle_endpoints(size_t start_index, size_t end_index);

        void reset_distances();

        void init_distances();

        void add_cycle_endpoint_pair(size_t start_index, size_t end_index);

        void add_cycle_endpoint_pair_mono(size_t start_index, size_t end_index);

        void dijkstra_shortest_distance_path_debug(size_t source_index, size_t target_index);

        void init_SO();

        void init_set_all_edge_weights();

        void reweight(pm::Mwpm &mwpm);

        void reweight_reset(pm::Mwpm &mwpm);

        void setting_edge_weights();

        int64_t SoftOutput(pm::Mwpm &mwpm);

        std::pair<int64_t, int64_t> SoftOutput_2d(pm::Mwpm &mwpm);

        int64_t dijkstra_shortest_distance_between(size_t source_index, size_t target_index);

        double normalizer;


        private:
        std::vector<DijkstraNodeVisitor> visitors;
        std::vector<vertex_descriptor> discovered_nodes;
        std::vector<int64_t> distances;
        std::vector<int64_t> cycle_lengths;
        int64_t inf;
        

    };


    void process_timeline_until_completion(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);
    size_t flooded_nodes_counter(pm::Mwpm& mwpm);
    size_t flooded_nodes_counter_by_region(pm::Mwpm& mwpm);

}

#endif
#include "dijkstra_graph.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include "pymatching/sparse_blossom/flooder/detector_node.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"

namespace dijkstra{
    
    VertexData::VertexData() : v_tag("inner") {}
    VertexData::VertexData(std::string v_tag) : v_tag(v_tag) {}


    WeightedEdgeData::WeightedEdgeData() : weight(0),
        weight_copy(0), flooded_weight_ts(nullptr), flooded_weight_st(nullptr) {}

    WeightedEdgeData::WeightedEdgeData(int64_t weight) : weight(weight),
        weight_copy(weight), flooded_weight_ts(nullptr), flooded_weight_st(nullptr) {} 


    DijkstraNodeVisitor::DijkstraNodeVisitor(const vertex_descriptor& target, 
                        std::vector<vertex_descriptor>& discovered_nodes) : 
                            target(target), discovered_nodes(discovered_nodes) {}


    int64_t get_final_local_radius_from_node(pm::DetectorNode *node) {
        if (node->region_that_arrived_top == nullptr) {
            return 0;
        }
        return node->region_that_arrived_top->radius.y_intercept() + node->wrapped_radius_cached;
    }
      

    SoftOutputDijkstra::SoftOutputDijkstra() {
        fl_matching_graph = flooded_graph_t();
        inf = std::numeric_limits<int64_t>::max();
        normalizer = inf;
    }

    SoftOutputDijkstra::SoftOutputDijkstra(pm::Mwpm &mwpm) {
        fl_matching_graph = mwpm_to_dijkstra_graph(mwpm);
        inf = std::numeric_limits<int64_t>::max();
        normalizer = mwpm.flooder.graph.normalising_constant;
    }

    void SoftOutputDijkstra::add_edge_from_mwpm_to_graph(size_t source_index, size_t target_index,
                                                pm::Mwpm &mwpm, 
                                                flooded_graph_t &graph) {
        int64_t edge_weight {0};
        int64_t *flooded_weight_ij {nullptr};
        int64_t *flooded_weight_ji {nullptr};
        


        pm::DetectorNode &source_node = mwpm.flooder.graph.nodes[source_index];
        for (size_t i = 0; i < source_node.neighbors.size(); i++) {
            if (source_node.neighbors[i]==nullptr) {
                continue;
            }
            if (target_index == (source_node.neighbors[i]-&mwpm.flooder.graph.nodes[0])) {
                edge_weight = nodes[source_index].neighbor_weights[i];
                if (source_index < target_index) {
                    flooded_weight_ij = &nodes[source_index].neighbor_flooded_weights[i];
                } else if (source_index > target_index)
                {
                    flooded_weight_ji = &nodes[source_index].neighbor_flooded_weights[i];
                }
                else {throw (2);}
                
            }
        }

        pm::DetectorNode &target_node = mwpm.flooder.graph.nodes[target_index];
        for (size_t i = 0; i < target_node.neighbors.size(); i++) {
            if (target_node.neighbors[i]==nullptr) {
                continue;
            }
            if (source_index == (target_node.neighbors[i]-&mwpm.flooder.graph.nodes[0])) {
                if (target_index < source_index) {
                    flooded_weight_ij = &nodes[target_index].neighbor_flooded_weights[i];
                } else if (target_index > source_index)
                {
                    flooded_weight_ji = &nodes[target_index].neighbor_flooded_weights[i];
                }
                else {throw (2);}
                
            }
        }

        vertex_descriptor source = boost::vertex(source_index, graph);
        vertex_descriptor target = boost::vertex(target_index, graph);
        edge_descriptor e;
        bool b;
        boost::tie(e,b) = boost::add_edge(source, target, WeightedEdgeData(edge_weight), graph);
        graph[e].flooded_weight_st = flooded_weight_ij;
        graph[e].flooded_weight_ts = flooded_weight_ji;
        // std::cout << *flooded_weight_ij << "==" << *flooded_weight_ji << "==" << graph[e].weight << std::endl;


    }

    flooded_graph_t SoftOutputDijkstra::mwpm_to_dijkstra_graph(pm::Mwpm &mwpm){
        for (size_t i = 0; i < mwpm.flooder.graph.num_nodes; i++) {
            Node node = Node();
            pm::DetectorNode &node_mwpm = mwpm.flooder.graph.nodes[i];
            for (size_t j = 0; j < node_mwpm.neighbor_weights.size(); j++) {
                node.neighbor_weights.push_back(node_mwpm.neighbor_weights[j]);
                node.neighbor_flooded_weights.push_back(node_mwpm.neighbor_weights[j]);
            }
            nodes.push_back(node);
        }
        flooded_graph_t weighted_graph = flooded_graph_t(mwpm.flooder.graph.num_nodes);
        for (pm::DetectorNode &node_source : mwpm.flooder.graph.nodes) {
            size_t source_index = &node_source - &mwpm.flooder.graph.nodes[0];
            vertex_descriptor source = boost::vertex(source_index, weighted_graph);
            for (size_t i = 0; i < node_source.neighbors.size() ; i++) {
                if (node_source.neighbors[i] == nullptr) {
                    continue;
                } 
                size_t target_index = node_source.neighbors[i] - &mwpm.flooder.graph.nodes[0];
                vertex_descriptor target = boost::vertex(target_index, weighted_graph);
                if (boost::edge(source, target, weighted_graph).second == false) {
                    add_edge_from_mwpm_to_graph(source_index, target_index, mwpm, weighted_graph);
                }
            }
        }
        return weighted_graph;
    }

    void SoftOutputDijkstra::init_distances() {
        distances.clear();
        distances.resize(boost::num_vertices(fl_matching_graph),inf);
        discovered_nodes.clear();
    }


    void SoftOutputDijkstra::reset_distances() {
        // pmap_v_index vertIndx = boost::get(boost::vertex_index,fl_matching_graph);
        // pmap_v_index vertIndx = boost::get(boost::vertex_index,fl_matching_graph);
        // for (vertex_descriptor &node : discovered_nodes) {
        //     distances[vertIndx[node]] = inf;
        // }
        std::fill(distances.begin(), distances.end(), inf);
        discovered_nodes.clear();
    }

    void SoftOutputDijkstra::add_boundary_node(size_t boundary_index) {
        if (boundary_index != boost::num_vertices(fl_matching_graph)) {
            std::cout << boost::num_vertices(fl_matching_graph);
            throw std::invalid_argument("boundary index incorrect; should be exactly 1 larger \
                                         than the current max vertex index in the graph");
        }
        boost::add_vertex(VertexData("boundary"),fl_matching_graph);
    }

    void SoftOutputDijkstra::add_boundary_edge_from_mwpm(size_t inner_index, size_t boundary_index,
                                                         pm::Mwpm &mwpm) {
        vertex_descriptor boundary_node = boost::vertex(boundary_index, fl_matching_graph);
        if (fl_matching_graph[boundary_node].v_tag != "boundary") {
            throw std::invalid_argument("incorrect boundary index. The corresponding node \
                                         is not a boundary node.");
        }
        vertex_descriptor inner_node = boost::vertex(inner_index, fl_matching_graph);
        pm::DetectorNode node_mwpm = mwpm.flooder.graph.nodes[inner_index];
        if (node_mwpm.neighbors.size() == 0) {
            throw std::invalid_argument("this node has no error mechanism..");
            return;}
        // if (inner_index == 300) {
        //     std::cout << "neighbor size" << node_mwpm.neighbors.size() << std::endl;
        //     std::cout << "neighbor index" << node_mwpm.neighbors[0] - &mwpm.flooder.graph.nodes[0] << std::endl;
        //     std::cout << "neighbor weight" << node_mwpm.neighbor_weights[0] << std::endl;
        // }
        if (node_mwpm.neighbors[0] != nullptr) {
            throw std::invalid_argument("incorrect inner index. No boundary edge is found \
                                         on this node.");
        }
        int64_t edge_weight = nodes[inner_index].neighbor_weights[0];
        edge_descriptor e; 
        bool b;
        boost::tie(e,b) = boost::add_edge(inner_node, boundary_node, WeightedEdgeData(edge_weight),fl_matching_graph);
        fl_matching_graph[e].flooded_weight_st = &nodes[inner_index].neighbor_flooded_weights[0];
        fl_matching_graph[e].flooded_weight_ts = fl_matching_graph[e].flooded_weight_st;
    }

    void SoftOutputDijkstra::add_boundary_edge_to_image_from_mwpm(size_t original_index, size_t image_index, size_t boundary_index, pm::Mwpm &mwpm) {
        vertex_descriptor boundary_node = boost::vertex(boundary_index, fl_matching_graph);
        if (fl_matching_graph[boundary_node].v_tag != "boundary") {
            throw std::invalid_argument("incorrect boundary index. The corresponding node \
                                         is not a boundary node.");
        }
        vertex_descriptor inner_node = boost::vertex(original_index, fl_matching_graph);
        vertex_descriptor image_node = boost::vertex(image_index, fl_matching_graph);
        pm::DetectorNode node_mwpm = mwpm.flooder.graph.nodes[original_index];
        if (node_mwpm.neighbors.size() == 0) {
            throw std::invalid_argument("this node has no error mechanism..");
            return;
        }
        // if (inner_index == 300) {
        //     std::cout << "neighbor size" << node_mwpm.neighbors.size() << std::endl;
        //     std::cout << "neighbor index" << node_mwpm.neighbors[0] - &mwpm.flooder.graph.nodes[0] << std::endl;
        //     std::cout << "neighbor weight" << node_mwpm.neighbor_weights[0] << std::endl;
        // }
        if (node_mwpm.neighbors[0] != nullptr) {
            throw std::invalid_argument("incorrect inner index. No boundary edge is found \
                                         on this node.");
        }
        int64_t edge_weight = nodes[original_index].neighbor_weights[0];
        edge_descriptor e; 
        bool b;
        boost::tie(e,b) = boost::add_edge(image_node, boundary_node, WeightedEdgeData(edge_weight),fl_matching_graph);
        fl_matching_graph[e].flooded_weight_st = &nodes[original_index].neighbor_flooded_weights[0];
        fl_matching_graph[e].flooded_weight_ts = fl_matching_graph[e].flooded_weight_st;
    }

    void SoftOutputDijkstra::copy_edge_to_image(size_t original_s_index, size_t original_t_index, size_t image_s_index, size_t image_t_index) {
        if (original_s_index >= boost::num_vertices(fl_matching_graph) || 
            original_t_index >= boost::num_vertices(fl_matching_graph) || 
            image_s_index >= boost::num_vertices(fl_matching_graph) || 
            image_t_index >= boost::num_vertices(fl_matching_graph)) {
            throw std::invalid_argument("index out of range");
        }
        vertex_descriptor original_s = boost::vertex(original_s_index,fl_matching_graph);
        vertex_descriptor image_s = boost::vertex(image_s_index,fl_matching_graph);
        vertex_descriptor original_t = boost::vertex(original_t_index,fl_matching_graph);
        vertex_descriptor image_t = boost::vertex(image_t_index,fl_matching_graph);
        edge_descriptor e;
        bool b;
        boost::tie(e,b) = boost::edge(original_s,original_t,fl_matching_graph);
        if (b==false) {
            throw std::invalid_argument("no edge between the original vertex and the target in the first place");
            return;
        }
        edge_descriptor e_image;
        bool b_image;
        boost::tie(e_image,b_image) = boost::add_edge(image_s,image_t,
                                                      WeightedEdgeData(fl_matching_graph[e].weight),
                                                      fl_matching_graph);
        if (b_image==false) {
            throw std::invalid_argument("failed attempt of redirecting the edge");
        }
        fl_matching_graph[e_image].flooded_weight_st = fl_matching_graph[e].flooded_weight_st;
        fl_matching_graph[e_image].flooded_weight_ts = fl_matching_graph[e].flooded_weight_ts;
    }

    void SoftOutputDijkstra::add_image_node(size_t original_index, size_t image_index) {
        if (original_index >= boost::num_vertices(fl_matching_graph)) {
            throw std::invalid_argument("inner index out of range");
        }
        if (image_index != boost::num_vertices(fl_matching_graph)) {
            throw std::invalid_argument("image index incorrect; should be exactly 1 larger" 
                                        "than the current max node index");
        }
        boost::add_vertex(VertexData("image"),fl_matching_graph);
    }

    void SoftOutputDijkstra::redirect_edge_to_image(size_t original_index, size_t image_index, size_t target_index) {
        if (original_index >= boost::num_vertices(fl_matching_graph) || 
            image_index >= boost::num_vertices(fl_matching_graph)) {
            throw std::invalid_argument("index out of range");
        }
        vertex_descriptor original = boost::vertex(original_index,fl_matching_graph);
        vertex_descriptor image = boost::vertex(image_index,fl_matching_graph);
        vertex_descriptor target = boost::vertex(target_index,fl_matching_graph);
        edge_descriptor e;
        bool b;
        boost::tie(e,b) = boost::edge(original,target,fl_matching_graph);
        if (b==false) {
            throw std::invalid_argument("no edge between the original vertex and the target in the first place");
            return;
        }
        edge_descriptor e_image;
        bool b_image;
        boost::tie(e_image,b_image) = boost::add_edge(image,target,
                                                      WeightedEdgeData(fl_matching_graph[e].weight),
                                                      fl_matching_graph);
        if (b_image==false) {
            throw std::invalid_argument("failed attempt of redirecting the edge");
        }
        fl_matching_graph[e_image].flooded_weight_st = fl_matching_graph[e].flooded_weight_st;
        fl_matching_graph[e_image].flooded_weight_ts = fl_matching_graph[e].flooded_weight_ts;
        boost::remove_edge(e,fl_matching_graph);
    }

    int64_t SoftOutputDijkstra::dijkstra_shortest_distance_between(size_t source_index, size_t target_index) {
        distances[source_index] = 0;
        DijkstraNodeVisitor vis = DijkstraNodeVisitor(boost::vertex(target_index,fl_matching_graph),
                                                        discovered_nodes);
        vertex_descriptor source = boost::vertex(source_index,fl_matching_graph);
        try {
            boost::dijkstra_shortest_paths_no_color_map(fl_matching_graph,source,
                            boost::distance_map(make_iterator_property_map(distances.begin(),
                                    boost::get(boost::vertex_index,fl_matching_graph)))
                            .weight_map(boost::get(&WeightedEdgeData::weight,fl_matching_graph))
                            .vertex_index_map(boost::get(boost::vertex_index,fl_matching_graph))
                            .visitor(vis).distance_inf(inf)
                            .distance_compare(std::less<int64_t>())
                            .distance_combine(boost::closed_plus<int64_t>()));
        }
        catch (int gotcha) {}
        // std::cout << "distance from " << source_index << " to " << target_index << ": ";
        // if (distances[target_index] < inf) {
        //     std::cout << distances[target_index]/normalizer << std::endl;
        // } 
        // else {
        //     std::cout << "not connected." << std::endl;
        // }
        int64_t ret_distance = distances[target_index];
        reset_distances();
        return ret_distance;       
    }

    void SoftOutputDijkstra::dijkstra_shortest_distance_path_debug(size_t source_index, size_t target_index) {
        init_distances();
        init_set_all_edge_weights();
        // std::cout << boost::num_vertices(fl_matching_graph) << std::endl;
        std::vector<vertex_descriptor> predecessor(boost::num_vertices(fl_matching_graph));
        vertex_iterator v, v_end;
        pmap_v_index vertIndx = boost::get(boost::vertex_index,fl_matching_graph);
        for (boost::tie(v,v_end) = boost::vertices(fl_matching_graph); v != v_end; v++) {
            predecessor[vertIndx(*v)] = *v;
        }
        distances[source_index] = 0;
        DijkstraNodeVisitor vis = DijkstraNodeVisitor(boost::vertex(target_index,fl_matching_graph),
                                                        discovered_nodes);
        vertex_descriptor source = boost::vertex(source_index,fl_matching_graph);
        try {
            boost::dijkstra_shortest_paths_no_color_map(fl_matching_graph,source,
                            boost::distance_map(make_iterator_property_map(distances.begin(),
                                    boost::get(boost::vertex_index,fl_matching_graph)))
                            .predecessor_map(&predecessor[0])
                            .weight_map(boost::get(&WeightedEdgeData::weight,fl_matching_graph))
                            .vertex_index_map(boost::get(boost::vertex_index,fl_matching_graph))
                            .visitor(vis).distance_inf(inf)
                            .distance_compare(std::less<int64_t>())
                            .distance_combine(boost::closed_plus<int64_t>()));
            // boost::dijkstra_shortest_paths_no_color_map_no_init(fl_matching_graph,source,
            //                 &predecessor[0],
            //                 boost::make_iterator_property_map(distances.begin(),
            //                                     boost::get(boost::vertex_index,fl_matching_graph)),
            //                 boost::get(&WeightedEdgeData::weight,fl_matching_graph),
            //                 boost::get(boost::vertex_index,fl_matching_graph),
            //                 std::less<int64_t>(),
            //                 boost::closed_plus<int64_t>(),
            //                 inf,0,vis);
        }
        catch (int gotcha) {}
        // // std::cout << "distance from " << source_index << " to " << target_index << ": ";
        // // if (distances[target_index] < inf) {
        // //     std::cout << distances[target_index]/normalizer << std::endl;
        // // } 
        // // else {
        // //     std::cout << "not connected." << std::endl;
        // // }
        int64_t ret_distance = distances[target_index];
        std::cout << ret_distance/normalizer << std::endl;
        // std::cout << 'distance from source: ' << source_index << ' to target ' << target_index 
        //           << ' is ' << ret_distance/normalizer << std::endl;
        
        // std::cout << 'target: ' << target_index;
        while (predecessor[target_index]!=boost::vertex(target_index,fl_matching_graph)) {
                std::cout << target_index << std::endl;
                target_index = vertIndx(predecessor[target_index]);
                std::cout << target_index << std::endl;
        }
        // std::cout << ' , source: ' << source_index << std::endl; 
        reset_distances();
    }

    
    void SoftOutputDijkstra::reweight(pm::Mwpm &mwpm) {
        flooded_node_indices.clear();
        // if (mwpm.flooder.region_arena.allocated.size() == 0) {
        //     return;
        // }
        // for (pm::GraphFillRegion *region : mwpm.flooder.region_arena.allocated) {
        //     for (pm::DetectorNode *node : region->shell_area) {
        //         size_t node_index = node - &mwpm.flooder.graph.nodes[0];
        //         if (node_index > nodes.size()-1) {
        //             std::cout << node_index << std::endl;
        //             throw std::invalid_argument("what? node index out of range.");
        //         }
                
        //         nodes[node_index].local_radius = get_final_local_radius_from_node(node);
        //         flooded_node_indices.push_back(node_index);
        //     }
        // }
        for (pm::DetectorNode &node : mwpm.flooder.graph.nodes) {
            if (node.region_that_arrived_top != nullptr) {
                size_t node_index = &node - &mwpm.flooder.graph.nodes[0];
                if (node_index > nodes.size()-1) {
                    std::cout << node_index << std::endl;
                    throw std::invalid_argument("what? node index out of range.");
                }
                nodes[node_index].local_radius = get_final_local_radius_from_node(&node);
                flooded_node_indices.push_back(node_index);
            }  
        }
        for (size_t node_index : flooded_node_indices) {
            pm::DetectorNode* node = &mwpm.flooder.graph.nodes[node_index];
            for (size_t i=0; i<node->neighbors.size(); i++) {
                if (node->neighbors[i] == nullptr) {
                    nodes[node_index].neighbor_flooded_weights[i] -= nodes[node_index].local_radius;
                    // node->neighbor_flooded_weights[i] -= node->local_radius_fixed; 
                    if (nodes[node_index].neighbor_flooded_weights[i] < 0) {
                        nodes[node_index].neighbor_flooded_weights[i] = 0;
                    }
                } else {
                    int64_t combined_radius {0};
                    combined_radius += nodes[node_index].local_radius;
                    size_t neighbor_index = node->neighbors[i] - &mwpm.flooder.graph.nodes[0];
                    combined_radius += nodes[neighbor_index].local_radius;
                    nodes[node_index].neighbor_flooded_weights[i] -= combined_radius; 
                    if (combined_radius < 0) {throw std::invalid_argument("combined_radius < 0");}
                    if (nodes[node_index].neighbor_flooded_weights[i] < 0) {
                        nodes[node_index].neighbor_flooded_weights[i] = 0;
                    }
                }
            }
        }
        setting_edge_weights();
    }

    int64_t SoftOutputDijkstra::SoftOutput(pm::Mwpm &mwpm){
        init_distances();
        reweight(mwpm);
        int64_t pathweight {inf};
        if (cycle_endpoint_pairs.size() == 0) {
            throw std::invalid_argument("cycle endpoints unset");
        }
        for (std::pair<size_t,size_t> cycle_endpoints : cycle_endpoint_pairs) {
            int64_t cycle_length = dijkstra_shortest_distance_between(cycle_endpoints.first,
                                                                      cycle_endpoints.second);
            if (pathweight > cycle_length) {
                pathweight = cycle_length;
            }
        }
        reweight_reset(mwpm);
        return pathweight;
    }

    std::pair<int64_t, int64_t> SoftOutputDijkstra::SoftOutput_2d(pm::Mwpm &mwpm){
        init_distances();
        reweight(mwpm);
        int64_t pathweight {inf};
        int64_t pathweight_mono {inf};
        if (cycle_endpoint_pairs.size() == 0) {
            throw std::invalid_argument("cycle endpoints unset");
        }
        for (std::pair<size_t,size_t> cycle_endpoints : cycle_endpoint_pairs) {
            int64_t cycle_length = dijkstra_shortest_distance_between(cycle_endpoints.first,
                                                                      cycle_endpoints.second);
            if (pathweight > cycle_length) {
                pathweight = cycle_length;
            }
        }
        if (cycle_endpoint_pairs_mono.size() == 0) {
            throw std::invalid_argument("cycle mono endpoints unset");
        }
        for (std::pair<size_t,size_t> cycle_endpoints : cycle_endpoint_pairs_mono) {
            int64_t cycle_length = dijkstra_shortest_distance_between(cycle_endpoints.first,
                                                                      cycle_endpoints.second);
            if (pathweight_mono > cycle_length) {
                pathweight_mono = cycle_length;
            }
        }
        reweight_reset(mwpm);
        return {pathweight_mono,pathweight};
    }

    void SoftOutputDijkstra::reweight_reset(pm::Mwpm &mwpm) {
        for (size_t node_index : flooded_node_indices) {
            nodes[node_index].local_radius = 0;
            for (size_t i=0; i<nodes[node_index].neighbor_weights.size(); i++) {
                nodes[node_index].neighbor_flooded_weights[i] = nodes[node_index].neighbor_weights[i];
            }
        }
        setting_edge_weights();
        flooded_node_indices.clear();
    }

    void SoftOutputDijkstra::setting_edge_weights() {
        for (size_t node_index : flooded_node_indices) {
            vertex_descriptor vertex = boost::vertex(node_index, fl_matching_graph);
            out_edge_iterator e_start, e_end;
            boost::tie(e_start, e_end) = boost::out_edges(vertex,fl_matching_graph);
            for (out_edge_iterator e_out = e_start ; e_out!=e_end; e_out++) {
                fl_matching_graph[*e_out].weight = std::min(*fl_matching_graph[*e_out].flooded_weight_st,
                                                            *fl_matching_graph[*e_out].flooded_weight_ts);
            }
        }
    }

    void SoftOutputDijkstra::init_set_all_edge_weights() {
        edge_iterator e, e_end;
        for (boost::tie(e,e_end)=boost::edges(fl_matching_graph); e!=e_end; e++) {
            fl_matching_graph[*e].weight = std::min(*fl_matching_graph[*e].flooded_weight_st,
                                                    *fl_matching_graph[*e].flooded_weight_ts);
        }
    }

    void SoftOutputDijkstra::add_cycle_endpoint_pair(size_t start_index, size_t end_index) {
        cycle_endpoint_pairs.push_back({start_index, end_index});
    }

    void SoftOutputDijkstra::add_cycle_endpoint_pair_mono(size_t start_index, size_t end_index) {
        cycle_endpoint_pairs_mono.push_back({start_index, end_index});
    }

    size_t flooded_nodes_counter(pm::Mwpm& mwpm) {
        size_t counter = 0;
        for (pm::DetectorNode& node : mwpm.flooder.graph.nodes) {
            if (node.region_that_arrived != nullptr) {
                counter++;
            }
        }
        return counter;
    }

    size_t flooded_nodes_counter_by_region(pm::Mwpm& mwpm) {
        size_t counter = 0;
        for (pm::GraphFillRegion *region : mwpm.flooder.region_arena.allocated) {
            counter += region->shell_area.size();
        }
        return counter;
    }

}
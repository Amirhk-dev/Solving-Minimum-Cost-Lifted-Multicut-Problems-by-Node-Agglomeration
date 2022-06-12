#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_CUTMIN_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_CUTMIN_HXX

#include <cstddef>
#include <iterator>
#include <vector>
#include <algorithm>
#include <map>
#include <queue>

#include "andres/partition.hxx"

namespace andres {
    namespace graph {
        namespace multicut_lifted {

            /// Greedy balanced agglomerative decomposition of a graph. (BEC)
            ///
            template<typename ORIGGRAPH, typename LIFTGRAPH, typename EVA, typename ELA>
                void balancedEdgeContraction(
                        const ORIGGRAPH& original_graph,
                        const LIFTGRAPH& lifted_graph,
                        EVA const& edge_values,
                        ELA& edge_labels
                        )
                { 
                    class DynamicGraph
                    {
                        public:
                            DynamicGraph(size_t n) :
                                vertices_(n),
                                vertex_weights_(n)
                        {}

                            bool edgeExists(size_t a, size_t b) const
                            {
                                return !vertices_[a].empty() && vertices_[a].find(b) != vertices_[a].end();
                            }

                            std::map<size_t, typename EVA::value_type> const& getAdjacentVertices(size_t v) const
                            {
                                return vertices_[v];
                            }

                            typename EVA::value_type getEdgeWeight(size_t a, size_t b) const
                            {
                                return vertices_[a].at(b);
                            }

                            void removeVertex(size_t v)
                            {
                                for (auto& p : vertices_[v])
                                    vertices_[p.first].erase(v);

                                vertices_[v].clear();
                            }

                            void setEdgeWeight(size_t a, size_t b, typename EVA::value_type w)
                            {
                                vertices_[a][b] = w;
                                vertices_[b][a] = w;
                            }

                            void setVertexWeights(size_t a, typename EVA::value_type w)
                            {
                                vertex_weights_[a] = w;  
                            }

                            typename EVA::value_type returnVertexWeights(size_t a)
                            {
                                return vertex_weights_[a];
                            }

                            void printVertexWeights(size_t a)
                            {
                                std::cout<< vertex_weights_[a]<<" "<<std::flush;
                            }

                        private:
                            std::vector<std::map<size_t, typename EVA::value_type>> vertices_;
                            std::vector<typename EVA::value_type> vertex_weights_;
                    };

                    struct Edge
                    {
                        Edge(size_t _a, size_t _b, typename EVA::value_type _wp, typename EVA::value_type _w)
                        {
                            if (_a > _b)
                                std::swap(_a, _b);

                            a = _a;
                            b = _b;
                            wp=_wp;
                            w = _w;
                        }

                        size_t a;
                        size_t b;
                        size_t edition;
                        typename EVA::value_type w;
                        typename EVA::value_type wp;

                        bool operator <(Edge const& other) const
                        {
                            return w < other.w;
                        }
                    };

                    std::vector<std::map<size_t, size_t>> edge_editions(original_graph.numberOfVertices());
                    DynamicGraph original_graph_cp(original_graph.numberOfVertices());   
                    DynamicGraph lifted_graph_cp(original_graph.numberOfVertices());

                    std::priority_queue<Edge> Q;

                    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
                    {
                        auto a = original_graph.vertexOfEdge(i, 0);
                        auto b = original_graph.vertexOfEdge(i, 1);

                        original_graph_cp.setEdgeWeight(a, b, 1.);
                    }


                    for (size_t i = 0; i < lifted_graph.numberOfVertices(); ++i)
                        lifted_graph_cp.setVertexWeights(i,1.0);


                    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
                    {
                        auto a = lifted_graph.vertexOfEdge(i, 0);
                        auto b = lifted_graph.vertexOfEdge(i, 1);
                        lifted_graph_cp.setEdgeWeight(a, b, edge_values[i]);
                        if (original_graph_cp.edgeExists(a, b))
                        {
                            auto e = Edge(a, b, edge_values[i],  edge_values[i]);
                            e.edition = ++edge_editions[e.a][e.b];
                            Q.push(e);
                        }
                    }

                    std::cout<<"graph initialized"<<std::endl;

                    andres::Partition<size_t> partition(original_graph.numberOfVertices());

                    float nmerges=0;
                    while (!Q.empty())
                    {
                        auto edge = Q.top();
                        Q.pop();

                        if (!original_graph_cp.edgeExists(edge.a, edge.b) || edge.edition < edge_editions[edge.a][edge.b])
                            continue;

                        if (edge.w < typename EVA::value_type())
                            break;

                        nmerges++;

                        auto stable_vertex = edge.a;
                        auto merge_vertex = edge.b;

                        if (lifted_graph_cp.getAdjacentVertices(stable_vertex).size() < lifted_graph_cp.getAdjacentVertices(merge_vertex).size())
                            std::swap(stable_vertex, merge_vertex);

                        partition.merge(stable_vertex, merge_vertex);

                        for (auto& p : original_graph_cp.getAdjacentVertices(merge_vertex))
                        {
                            if (p.first == stable_vertex)
                                continue;

                            original_graph_cp.setEdgeWeight(stable_vertex, p.first, 1.);
                        }

                        original_graph_cp.removeVertex(merge_vertex);
                        auto nwa = typename EVA::value_type(); 
                        auto nwb = typename EVA::value_type();
                        nwa = lifted_graph_cp.returnVertexWeights(stable_vertex);
                        nwb = lifted_graph_cp.returnVertexWeights(merge_vertex);

                        lifted_graph_cp.setVertexWeights(stable_vertex,nwa + nwb);

                        for (auto& p : lifted_graph_cp.getAdjacentVertices(stable_vertex))
                        {
                            if (p.first == merge_vertex)
                                continue;

                            auto nwp = typename EVA::value_type();
                            nwp =  lifted_graph_cp.returnVertexWeights(p.first);

                            auto tp = typename EVA::value_type();
                            if(lifted_graph_cp.edgeExists(merge_vertex, p.first))
                                continue;

                            float wn = (nwa+nwb+nwp)/(static_cast<float>(original_graph.numberOfVertices())/nmerges);
                            if (original_graph_cp.edgeExists(stable_vertex, p.first))
                            {
                                auto e = Edge(stable_vertex, p.first, p.second + tp, (p.second +tp)/wn);

                                e.edition = ++edge_editions[e.a][e.b];

                                Q.push(e);
                            }
                        }

                        for (auto& p : lifted_graph_cp.getAdjacentVertices(merge_vertex))
                        {
                            if (p.first == stable_vertex)
                                continue;

                            auto nwp = typename EVA::value_type();
                            nwp =  lifted_graph_cp.returnVertexWeights(p.first);

                            auto tp = typename EVA::value_type();
                            if(lifted_graph_cp.edgeExists(stable_vertex, p.first))
                                tp = lifted_graph_cp.getEdgeWeight(stable_vertex, p.first);
                            else
                                tp=0;
                            lifted_graph_cp.setEdgeWeight(stable_vertex, p.first, p.second + tp);

                            float wn = (nwa+nwb+nwp)/(static_cast<float>(original_graph.numberOfVertices())/nmerges);
                            if (original_graph_cp.edgeExists(stable_vertex, p.first))
                            {
                                auto e = Edge(stable_vertex, p.first, p.second + tp, (p.second +tp)/wn);

                                e.edition = ++edge_editions[e.a][e.b];

                                Q.push(e);
                            }
                        }

                        lifted_graph_cp.removeVertex(merge_vertex);
                    }

                    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
                        edge_labels[i] = partition.find(lifted_graph.vertexOfEdge(i, 0)) == partition.find(lifted_graph.vertexOfEdge(i, 1)) ? 0 : 1;
                }

        } // namespace multicut_lifted 
    } // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_CUTMIN_HXX
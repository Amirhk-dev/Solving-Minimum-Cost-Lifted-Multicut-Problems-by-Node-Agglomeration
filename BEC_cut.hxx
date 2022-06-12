#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_MINMAX_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_MINMAX_HXX

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

            /// Greedy agglomerative balanced min-max decomposition of a graph. (BEC-cut)
            ///
            template<typename ORIGGRAPH, typename LIFTGRAPH, typename EVA, typename ELA>
                void balancedEdgeContraction_cut(
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
                            if(w==other.w)
                                return wp< other.wp;
                            else 
                                return w < other.w;
                        }
                    };

                    std::vector<std::map<size_t, size_t>> edge_editions(original_graph.numberOfVertices());
                    DynamicGraph original_graph_cp(original_graph.numberOfVertices());
                    DynamicGraph dual_graph_cp(original_graph.numberOfVertices());
                    DynamicGraph lifted_graph_cp(original_graph.numberOfVertices());

                    std::priority_queue<Edge> Q;

                    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
                    {
                        auto a = original_graph.vertexOfEdge(i, 0);
                        auto b = original_graph.vertexOfEdge(i, 1);

                        original_graph_cp.setEdgeWeight(a, b, 1.);
                    }

                    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
                    {
                        auto a = lifted_graph.vertexOfEdge(i, 0);
                        auto b = lifted_graph.vertexOfEdge(i, 1);

                        lifted_graph_cp.setEdgeWeight(a, b, edge_values[i]);
                    }

                    for (size_t i = 0; i < lifted_graph.numberOfVertices(); ++i)
                        dual_graph_cp.setVertexWeights(i,0);

                    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
                    {
                        auto a = lifted_graph.vertexOfEdge(i, 0);
                        auto b = lifted_graph.vertexOfEdge(i, 1);
                        auto dual_nw = typename EVA::value_type();

                        dual_nw = dual_graph_cp.returnVertexWeights(a);	
                        dual_graph_cp.setVertexWeights(a,dual_nw+edge_values[i]);
                        dual_nw = dual_graph_cp.returnVertexWeights(b);
                        dual_graph_cp.setVertexWeights(b,dual_nw+edge_values[i]);
                    }

                    for (size_t i = 0; i < lifted_graph.numberOfVertices(); ++i)
                        lifted_graph_cp.setVertexWeights(i,1.0);

                    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
                    {
                        auto a = lifted_graph.vertexOfEdge(i, 0);
                        auto b = lifted_graph.vertexOfEdge(i, 1);
                        auto dual_ew = typename EVA::value_type();
                        auto dual_nwa = typename EVA::value_type();
                        auto dual_nwb = typename EVA::value_type();
                        dual_nwa = dual_graph_cp.returnVertexWeights(a);
                        dual_nwb = dual_graph_cp.returnVertexWeights(b);

                        dual_ew= -(dual_nwa + dual_nwb - 2*edge_values[i]);

                        if (original_graph_cp.edgeExists(a, b))
                        {
                            auto e = Edge(a, b, dual_ew,  edge_values[i]);
                            e.edition = ++edge_editions[e.a][e.b];
                            Q.push(e);

                            if(edge_values[i]>=0){
                                e = Edge(a, b,   edge_values[i],dual_ew);
                                e.edition = edge_editions[e.a][e.b];
                            }
                        }
                    }

                    std::cout<<"graph initialized"<<std::endl;

                    andres::Partition<size_t> partition(original_graph.numberOfVertices());
                    float nmerges=0;

                    while (!(Q.empty()) )
                    {
                        auto edge = Q.top();
                        Q.pop();
                        if (edge.edition < edge_editions[edge.a][edge.b])
                            continue;

                        if (!original_graph_cp.edgeExists(edge.a, edge.b))
                            continue;

                        if (edge.w< typename EVA::value_type()&&edge.wp < typename EVA::value_type())
                        {
                            break;
                        }

                        if(edge.w < typename EVA::value_type())
                        {    
                            break;
                        }

                        nmerges++;

                        auto stable_vertex = edge.a;
                        auto merge_vertex = edge.b;

                        if (lifted_graph_cp.getAdjacentVertices(stable_vertex).size()
                                < lifted_graph_cp.getAdjacentVertices(merge_vertex).size())
                            std::swap(stable_vertex, merge_vertex);

                        partition.merge(stable_vertex, merge_vertex);

                        for (auto& p : original_graph_cp.getAdjacentVertices(merge_vertex))
                        {
                            if (p.first == stable_vertex)
                                continue;

                            original_graph_cp.setEdgeWeight(stable_vertex, p.first, 1.);
                        }

                        original_graph_cp.removeVertex(merge_vertex);

                        auto dual_nwa = dual_graph_cp.returnVertexWeights(stable_vertex);
                        auto dual_nwb = dual_graph_cp.returnVertexWeights(merge_vertex);
                        dual_nwa= dual_nwa + dual_nwb - 2*edge.wp;
                        dual_graph_cp.setVertexWeights(stable_vertex,dual_nwa);

                        auto nwa = lifted_graph_cp.returnVertexWeights(stable_vertex);
                        auto nwb = lifted_graph_cp.returnVertexWeights(merge_vertex);
                        lifted_graph_cp.setVertexWeights(stable_vertex,nwa + nwb);

                        for (auto& p : lifted_graph_cp.getAdjacentVertices(stable_vertex))
                        { 
                            if (p.first == merge_vertex)
                                continue;

                            if (!original_graph_cp.edgeExists(stable_vertex, p.first))
                                continue;

                            if(lifted_graph_cp.edgeExists(merge_vertex, p.first))
                                continue;

                            auto dual_nwp = dual_graph_cp.returnVertexWeights(p.first);
                            auto nwp =  lifted_graph_cp.returnVertexWeights(p.first);

                            float wn = (nwa+nwb+nwp)/(static_cast<float>(original_graph.numberOfVertices())/nmerges);

                            auto t=-(dual_nwa + dual_nwp - 2*p.second)/wn;

                            auto e = Edge(stable_vertex, p.first, t, p.second/wn  );
                            e.edition = ++edge_editions[e.a][e.b];
                            Q.push(e);

                            if(p.second<0)
                                continue;

                            e = Edge(stable_vertex, p.first, p.second/wn, t);
                            e.edition = edge_editions[e.a][e.b];
                        }

                        for (auto& p : lifted_graph_cp.getAdjacentVertices(merge_vertex))
                        {
                            if (p.first == stable_vertex)
                                continue;

                            auto dual_nwp = dual_graph_cp.returnVertexWeights(p.first);
                            auto nwp =  lifted_graph_cp.returnVertexWeights(p.first);

                            auto tp = typename EVA::value_type();
                            if(lifted_graph_cp.edgeExists(stable_vertex, p.first))
                                tp = lifted_graph_cp.getEdgeWeight(stable_vertex, p.first);
                            else
                                tp=0;

                            lifted_graph_cp.setEdgeWeight(stable_vertex, p.first, p.second + tp);
                            float wn = (nwa+nwb+nwp)/(static_cast<float>(original_graph.numberOfVertices())/nmerges);

                            auto t=-(dual_nwa + dual_nwp - 2*tp - 2*p.second)/wn;

                            if (original_graph_cp.edgeExists(stable_vertex, p.first))
                            {
                                auto e = Edge(stable_vertex, p.first, t,(p.second +tp)/wn );
                                e.edition = ++edge_editions[e.a][e.b];
                                Q.push(e);
                                if(p.second+tp<0)
                                    continue;
                                e = Edge(stable_vertex, p.first, (p.second +tp)/wn, t);
                                e.edition = edge_editions[e.a][e.b];
                            }
                        }

                        dual_graph_cp.removeVertex(merge_vertex);
                        lifted_graph_cp.removeVertex(merge_vertex);
                    }

                    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
                        edge_labels[i] = partition.find(lifted_graph.vertexOfEdge(i, 0)) == partition.find(lifted_graph.vertexOfEdge(i, 1)) ? 0 : 1;
                }

        } // namespace multicut_lifted 
    } // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_MINMAX_HXX
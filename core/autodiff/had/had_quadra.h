/**
    HAD is a single header C++ reverse-mode automatic differentiation library using operator overloading, with focus on
    second-order derivatives (Hessian).

    Quadra extension note:
    This header is an extension of the original had.h code by Tzu-Mao Li. It preserves the original MIT-licensed
    reverse-mode and edge-pushing Hessian implementation, and has been extended and modified to improve the
    Quadra framework by adding experimental third-order derivative support. The current third-order API focuses on
    exact directional third derivatives through a lightweight third-order forward scalar, while retaining the original
    reverse-mode Hessian machinery for value, gradient, and Hessian calculations.
    It implements the edge_pushing algorithm (see "Hessian Matrices via Automatic Differentiation",
    Gower and Mello 2010) to efficiently compute the second derivatives.

    See https://github.com/BachiLi/had for more details.

    Author: Tzu-Mao Li


    The MIT License (MIT)

    Copyright (c) 2015 Tzu-Mao Li

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
**/

#ifndef HAD_QUADRA_H__
#define HAD_QUADRA_H__

#include <vector>
#include <cmath>
#ifdef WIN32
#define threadDefine thread_local
#endif
#ifdef __APPLE__
#define USE_AATREE
#define threadDefine __thread
#endif

#ifdef __unix
// #define USE_AATREE
#define threadDefine __thread
#endif

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <utility>

#ifndef M_PI
#define M_PI std::acos(-1)
#endif

namespace had
{

    // Change the following line if you want to use single precision floats
    typedef double Real;
    typedef unsigned int VertexId;

    struct ADGraph;
    struct AReal;

    extern threadDefine ADGraph *g_ADGraph;
// Declare this in your .cpp source
#define DECLARE_ADGRAPH()                    \
    namespace had                            \
    {                                        \
        threadDefine ADGraph *g_ADGraph = 0; \
    }

    AReal NewAReal(const Real val);

    struct AReal
    {
        AReal() {}

        AReal(const Real val)
        {
            *this = NewAReal(val);
        }

        AReal(const Real val, const VertexId varId) : val(val), dot(Real(0.0)), varId(varId) {}

        Real val;
        // First-order directional tangent used by directional edge-pushing.
        // Important: this must be propagated on the AReal object itself,
        // not only in ADGraph::vertices, because later operator overloads
        // read operands' .dot values directly.
        Real dot = Real(0.0);
        VertexId varId;
    };

    struct ADEdge
    {
        ADEdge() {}
        ADEdge(const VertexId to, const Real w = Real(0.0), const Real dw = Real(0.0))
            : to(to), w(w), dw(dw) {}

        VertexId to;
        Real w;
        Real dw; // directional derivative of edge weight
    };

    // We assume there is at most 2 outgoing edges from this vertex
    struct ADVertex
    {
        ADVertex(const VertexId newId)
        {
            e1 = e2 = ADEdge(newId);
            w = wDot = soW = soWDot = toW = dot = Real(0.0);
        }

        // If ei.to == the id of this vertex, then the edge does not exist
        ADEdge e1, e2;
        // first-order adjoint weight
        Real w;
        // directional derivative of first-order adjoint weight
        Real wDot;
        // second-order weights
        // for vertex with single outgoing edge,
        // soW represents the second-order weight of the conntecting vertex (d^2f/dx^2)
        // for vertex with two outgoing edges,
        // soW represents the second-order weight between the conntecting vertices (d^2f/dxdy)
        // the system assumes d^2f/dx^2 & d^2f/dy^2 are both zero in the two outgoing edges case to save memory
        Real soW;
        // directional derivative of soW along seeded primal tangent
        Real soWDot;
        // third-order local derivative weight. For unary vertices this is d^3 child / d parent^3.
        // For binary vertices this is reserved for future full third-order edge-pushing support.
        Real toW;
        // Optional directional tangent associated with this vertex.
        Real dot;
    };

    struct BTNode
    {
        BTNode() {}
        BTNode(const VertexId key, const Real val) : key(key), val(val)
        {
            left = right = -1;
#ifdef USE_AATREE
            level = 1;
#endif
        }

        VertexId key;
        Real val;
        int left;
        int right;
#ifdef USE_AATREE
        int level;
#endif
    };

    struct BTree
    {
        BTree()
        {
            nodes.reserve(32);
            root = 0;
        }
#ifdef USE_AATREE
        inline void Skew()
        {
            if (nodes.size() == 0)
                return;

            while (nodes[root].left != -1 &&
                   nodes[nodes[root].left].level == nodes[root].level)
            {
                int l = nodes[root].left;
                nodes[root].left = nodes[l].right;
                nodes[l].right = root;
                root = l;
            }
        }

        inline void Split()
        {
            if (nodes.size() == 0)
                return;

            while (nodes[root].right != -1 &&
                   nodes[nodes[root].right].right != -1 &&
                   nodes[root].level == nodes[nodes[nodes[root].right].right].level)
            {
                int r = nodes[root].right;
                nodes[root].right = nodes[r].left;
                nodes[r].left = root;
                nodes[r].level++;
                root = r;
            }
        }
#endif
        inline void Insert(const VertexId key, const Real val)
        {
            int index = root;
            if (nodes.size() > 0)
            {
                int *lastEdge;
                do
                {
                    if (key == nodes[index].key)
                    {
                        nodes[index].val += val;
                        return;
                    }
                    lastEdge = &(nodes[index].left) + (key > nodes[index].key);
                    index = *lastEdge;
                } while (index >= 0);

                *lastEdge = nodes.size();
            }
            nodes.push_back(BTNode(key, val));
#ifdef USE_AATREE
            Skew();
            Split();
#endif
        }

        inline Real Query(const VertexId key)
        {
            int index = root;
            while (index >= 0 && index < (int)nodes.size())
            {
                if (key == nodes[index].key)
                {
                    return nodes[index].val;
                }
                else if (key < nodes[index].key)
                {
                    index = nodes[index].left;
                }
                else
                {
                    index = nodes[index].right;
                }
            }
            return Real(0.0);
        }

        inline void Clear()
        {
            nodes.clear();
            root = 0;
        }

        std::vector<BTNode> nodes;
        int root;
    };

    struct ADGraph
    {
        ADGraph()
        {
            g_ADGraph = this;
        }

        inline void Clear()
        {
            vertices.clear();
            soEdges.clear();
            selfSoEdges.clear();
            soEdgesDot.clear();
            selfSoEdgesDot.clear();
        }

        std::vector<ADVertex> vertices;
        std::vector<BTree> soEdges;
        std::vector<Real> selfSoEdges;
        std::vector<BTree> soEdgesDot;
        std::vector<Real> selfSoEdgesDot;
    };

    inline AReal NewAReal(const Real val)
    {
        std::vector<ADVertex> &vertices = g_ADGraph->vertices;
        VertexId newId = vertices.size();
        vertices.push_back(ADVertex(newId));
        return AReal(val, newId);
    }

    inline void AddEdge(AReal &c, const AReal &p,
                        const Real w, const Real soW, const Real toW = Real(0.0))
    {
        ADVertex &v = g_ADGraph->vertices[c.varId];
        const Real dw = soW * p.dot;
        v.e1 = ADEdge(p.varId, w, dw);
        v.soW = soW;
        v.soWDot = toW * p.dot;
        v.toW = toW;
        v.dot = w * p.dot;
        c.dot = v.dot;
    }
    inline void AddEdge(AReal &c,
                        const AReal &p1, const AReal &p2,
                        const Real w1, const Real w2,
                        const Real soW, const Real toW = Real(0.0))
    {
        ADVertex &v = g_ADGraph->vertices[c.varId];

        Real dw1 = Real(0.0);
        Real dw2 = Real(0.0);
        if (soW != Real(0.0))
        {
            dw1 = soW * p2.dot;
            dw2 = soW * p1.dot;
        }

        v.e1 = ADEdge(p1.varId, w1, dw1);
        v.e2 = ADEdge(p2.varId, w2, dw2);
        v.soW = soW;
        v.soWDot = toW * (p1.dot + p2.dot);
        v.toW = toW;
        v.dot = w1 * p1.dot + w2 * p2.dot;
        c.dot = v.dot;
    }

    ////////////////////// Addition ///////////////////////////
    inline AReal operator+(const AReal &l, const AReal &r)
    {
        AReal ret = NewAReal(l.val + r.val);
        AddEdge(ret, l, r, Real(1.0), Real(1.0), Real(0.0));
        return ret;
    }
    inline AReal operator+(const AReal &l, const Real r)
    {
        AReal ret = NewAReal(l.val + r);
        AddEdge(ret, l, Real(1.0), Real(0.0));
        return ret;
    }
    inline AReal operator+(const Real l, const AReal &r)
    {
        return r + l;
    }
    inline AReal &operator+=(AReal &l, const AReal &r)
    {
        return (l = l + r);
    }
    inline AReal &operator+=(AReal &l, const Real r)
    {
        return (l = l + r);
    }
    ///////////////////////////////////////////////////////////

    ////////////////// Subtraction ////////////////////////////
    inline AReal operator-(const AReal &l, const AReal &r)
    {
        AReal ret = NewAReal(l.val - r.val);
        AddEdge(ret, l, r, Real(1.0), -Real(1.0), Real(0.0));
        return ret;
    }
    inline AReal operator-(const AReal &l, const Real r)
    {
        AReal ret = NewAReal(l.val - r);
        AddEdge(ret, l, Real(1.0), Real(0.0));
        return ret;
    }
    inline AReal operator-(const Real l, const AReal &r)
    {
        AReal ret = NewAReal(l - r.val);
        AddEdge(ret, r, Real(-1.0), Real(0.0));
        return ret;
    }
    inline AReal &operator-=(AReal &l, const AReal &r)
    {
        return (l = l - r);
    }
    inline AReal &operator-=(AReal &l, const Real r)
    {
        return (l = l - r);
    }
    inline AReal operator-(const AReal &x)
    {
        AReal ret = NewAReal(-x.val);
        AddEdge(ret, x, Real(-1.0), Real(0.0));
        return ret;
    }
    ///////////////////////////////////////////////////////////

    ////////////////// Multiplication /////////////////////////
    inline AReal operator*(const AReal &l, const AReal &r)
    {
        AReal ret = NewAReal(l.val * r.val);
        AddEdge(ret, l, r, r.val, l.val, Real(1.0));
        return ret;
    }
    inline AReal operator*(const AReal &l, const Real r)
    {
        AReal ret = NewAReal(l.val * r);
        AddEdge(ret, l, r, Real(0.0));
        return ret;
    }
    inline AReal operator*(const Real l, const AReal &r)
    {
        return r * l;
    }
    inline AReal &operator*=(AReal &l, const AReal &r)
    {
        return (l = l * r);
    }
    inline AReal &operator*=(AReal &l, const Real r)
    {
        return (l = l * r);
    }
    ///////////////////////////////////////////////////////////

    ////////////////// Inversion //////////////////////////////
    inline AReal Inv(const AReal &x)
    {
        Real invX = Real(1.0) / x.val;
        Real invXSq = invX * invX;
        Real invXCu = invXSq * invX;
        AReal ret = NewAReal(invX);
        AddEdge(ret, x, -invXSq, Real(2.0) * invXCu, -Real(6.0) * invXCu * invX);
        return ret;
    }
    inline Real Inv(const Real x)
    {
        return Real(1.0) / x;
    }
    ///////////////////////////////////////////////////////////

    ////////////////// Division ///////////////////////////////
    inline AReal operator/(const AReal &l, const AReal &r)
    {
        return l * Inv(r);
    }
    inline AReal operator/(const AReal &l, const Real r)
    {
        return l * Inv(r);
    }
    inline AReal operator/(const Real l, const AReal &r)
    {
        return l * Inv(r);
    }
    inline AReal &operator/=(AReal &l, const AReal &r)
    {
        return (l = l / r);
    }
    inline AReal &operator/=(AReal &l, const Real r)
    {
        return (l = l / r);
    }
    ///////////////////////////////////////////////////////////

    ////////////////// Comparisons ////////////////////////////
    inline bool operator<(const AReal &l, const AReal &r)
    {
        return l.val < r.val;
    }
    inline bool operator<=(const AReal &l, const AReal &r)
    {
        return l.val <= r.val;
    }
    inline bool operator>(const AReal &l, const AReal &r)
    {
        return l.val > r.val;
    }
    inline bool operator>=(const AReal &l, const AReal &r)
    {
        return l.val >= r.val;
    }
    inline bool operator==(const AReal &l, const AReal &r)
    {
        return l.val == r.val;
    }
    ///////////////////////////////////////////////////////////

    //////////////// Misc functions ///////////////////////////
    inline Real square(const Real x)
    {
        return x * x;
    }
    inline AReal square(const AReal &x)
    {
        Real sqX = x.val * x.val;
        AReal ret = NewAReal(sqX);
        AddEdge(ret, x, Real(2.0) * x.val, Real(2.0), Real(0.0));
        return ret;
    }
    inline AReal sqrt(const AReal &x)
    {
        Real sqrtX = std::sqrt(x.val);
        Real invSqrtX = Real(1.0) / sqrtX;
        AReal ret = NewAReal(sqrtX);
        AddEdge(ret, x, Real(0.5) * invSqrtX, -Real(0.25) * invSqrtX / x.val, Real(0.375) * invSqrtX / (x.val * x.val));
        return ret;
    }
    inline AReal pow(const AReal &x, const Real a)
    {
        Real powX = std::pow(x.val, a);
        AReal ret = NewAReal(powX);
        AddEdge(ret, x, a * std::pow(x.val, a - Real(1.0)),
                a * (a - Real(1.0)) * std::pow(x.val, a - Real(2.0)),
                a * (a - Real(1.0)) * (a - Real(2.0)) * std::pow(x.val, a - Real(3.0)));
        return ret;
    }
    inline AReal exp(const AReal &x)
    {
        Real expX = std::exp(x.val);
        AReal ret = NewAReal(expX);
        AddEdge(ret, x, expX, expX, expX);
        return ret;
    }
    inline AReal log(const AReal &x)
    {
        Real logX = std::log(x.val);
        AReal ret = NewAReal(logX);
        Real invX = Real(1.0) / x.val;
        AddEdge(ret, x, invX, -invX * invX, Real(2.0) * invX * invX * invX);
        return ret;
    }
    inline AReal sin(const AReal &x)
    {
        Real sinX = std::sin(x.val);
        AReal ret = NewAReal(sinX);
        AddEdge(ret, x, std::cos(x.val), -sinX, -std::cos(x.val));
        return ret;
    }
    inline AReal cos(const AReal &x)
    {
        Real cosX = std::cos(x.val);
        AReal ret = NewAReal(cosX);
        AddEdge(ret, x, -std::sin(x.val), -cosX, std::sin(x.val));
        return ret;
    }
    inline AReal tan(const AReal &x)
    {
        Real tanX = std::tan(x.val);
        Real secX = Real(1.0) / std::cos(x.val);
        Real sec2X = secX * secX;
        AReal ret = NewAReal(tanX);
        AddEdge(ret, x, sec2X, Real(2.0) * tanX * sec2X,
                Real(2.0) * sec2X * sec2X + Real(4.0) * tanX * tanX * sec2X);
        return ret;
    }
    inline AReal asin(const AReal &x)
    {
        Real asinX = std::asin(x.val);
        AReal ret = NewAReal(asinX);
        Real tmp = Real(1.0) / (Real(1.0) - x.val * x.val);
        Real sqrtTmp = std::sqrt(tmp);
        AddEdge(ret, x, sqrtTmp, x.val * sqrtTmp * tmp,
                (Real(1.0) + Real(2.0) * x.val * x.val) * sqrtTmp * tmp * tmp);
        return ret;
    }
    inline AReal acos(const AReal &x)
    {
        Real acosX = std::acos(x.val);
        AReal ret = NewAReal(acosX);
        Real tmp = Real(1.0) / (Real(1.0) - x.val * x.val);
        Real negSqrtTmp = -std::sqrt(tmp);
        AddEdge(ret, x, negSqrtTmp, x.val * negSqrtTmp * tmp,
                -(Real(1.0) + Real(2.0) * x.val * x.val) * std::sqrt(tmp) * tmp * tmp);
        return ret;
    }
    ///////////////////////////////////////////////////////////

    inline void SetAdjoint(const AReal &v, const Real adj)
    {
        g_ADGraph->vertices[v.varId].w = adj;
    }

    inline Real GetAdjoint(const AReal &v)
    {
        return g_ADGraph->vertices[v.varId].w;
    }

    inline Real GetAdjoint(const AReal &i, const AReal &j)
    {
        if (i.varId == j.varId)
        {
            return g_ADGraph->selfSoEdges[i.varId];
        }
        else
        {
            return g_ADGraph->soEdges[std::max(i.varId, j.varId)].Query(std::min(i.varId, j.varId));
        }
    }

    inline VertexId SingleEdgePropagate(VertexId x, Real &a)
    {
        bool cont = g_ADGraph->vertices[x].e1.to != x &&
                    g_ADGraph->vertices[x].e2.to == x;
        while (cont)
        {
            a *= g_ADGraph->vertices[x].e1.w;
            x = g_ADGraph->vertices[x].e1.to;
            cont = g_ADGraph->vertices[x].e1.to != x &&
                   g_ADGraph->vertices[x].e2.to == x;
        }
        return x;
    }

    inline void PushEdge(const ADEdge &foEdge, const ADEdge &soEdge)
    {
        if (foEdge.to == soEdge.to)
        {
            g_ADGraph->selfSoEdges[foEdge.to] += Real(2.0) * foEdge.w * soEdge.w;
        }
        else
        {
            g_ADGraph->soEdges[std::max(foEdge.to, soEdge.to)].Insert(
                std::min(foEdge.to, soEdge.to), foEdge.w * soEdge.w);
        }
    }


    inline void PushEdgeDot(const ADEdge &foEdge,
                            const ADEdge &soEdge,
                            const Real soEdgeDot)
    {
        const Real valDot =
            foEdge.dw * soEdge.w +
            foEdge.w * soEdgeDot;

        if (foEdge.to == soEdge.to)
        {
            g_ADGraph->selfSoEdgesDot[foEdge.to] += Real(2.0) * valDot;
        }
        else
        {
            g_ADGraph->soEdgesDot[std::max(foEdge.to, soEdge.to)].Insert(
                std::min(foEdge.to, soEdge.to),
                valDot);
        }
    }

    inline Real GetAdjointDot(const AReal &i, const AReal &j)
    {
        if (i.varId == j.varId)
        {
            return g_ADGraph->selfSoEdgesDot[i.varId];
        }
        return g_ADGraph->soEdgesDot[std::max(i.varId, j.varId)].Query(
            std::min(i.varId, j.varId));
    }

    inline void PropagateAdjointDirectional()
    {
        const VertexId n_vertices =
            static_cast<VertexId>(g_ADGraph->vertices.size());

        if (g_ADGraph->soEdges.size() < g_ADGraph->vertices.size())
        {
            g_ADGraph->soEdges.resize(g_ADGraph->vertices.size());
        }
        else
        {
            for (int i = 0; i < (int)g_ADGraph->soEdges.size(); ++i)
                g_ADGraph->soEdges[i].Clear();
        }

        if (g_ADGraph->soEdgesDot.size() < g_ADGraph->vertices.size())
        {
            g_ADGraph->soEdgesDot.resize(g_ADGraph->vertices.size());
        }
        else
        {
            for (int i = 0; i < (int)g_ADGraph->soEdgesDot.size(); ++i)
                g_ADGraph->soEdgesDot[i].Clear();
        }

        g_ADGraph->selfSoEdges.assign(g_ADGraph->vertices.size(), Real(0.0));
        g_ADGraph->selfSoEdgesDot.assign(g_ADGraph->vertices.size(), Real(0.0));

        for (VertexId vid = n_vertices - 1; vid > 0; --vid)
        {
            ADVertex &vertex = g_ADGraph->vertices[vid];
            ADEdge &e1 = vertex.e1;
            ADEdge &e2 = vertex.e2;

            if (e1.to == vid)
                continue;

            // Push sparse off-diagonal Hessian edges and directional edges.
            BTree &btree = g_ADGraph->soEdges[vid];
            BTree &btreeDot = g_ADGraph->soEdgesDot[vid];

            for (auto it = btree.nodes.begin(); it != btree.nodes.end(); ++it)
            {
                ADEdge soEdge(it->key, it->val);
                const Real soDot = btreeDot.Query(it->key);

                PushEdge(e1, soEdge);
                PushEdgeDot(e1, soEdge, soDot);

                if (e2.to != vid)
                {
                    PushEdge(e2, soEdge);
                    PushEdgeDot(e2, soEdge, soDot);
                }
            }

            // Push diagonal Hessian entry.
            const Real S = g_ADGraph->selfSoEdges[vid];
            const Real SDot = g_ADGraph->selfSoEdgesDot[vid];

            if (S != Real(0.0) || SDot != Real(0.0))
            {
                g_ADGraph->selfSoEdges[e1.to] += e1.w * e1.w * S;
                g_ADGraph->selfSoEdgesDot[e1.to] +=
                    Real(2.0) * e1.w * e1.dw * S +
                    e1.w * e1.w * SDot;

                if (e2.to != vid)
                {
                    g_ADGraph->selfSoEdges[e2.to] += e2.w * e2.w * S;
                    g_ADGraph->selfSoEdgesDot[e2.to] +=
                        Real(2.0) * e2.w * e2.dw * S +
                        e2.w * e2.w * SDot;

                    const Real cross = e1.w * e2.w * S;
                    const Real crossDot =
                        (e1.dw * e2.w + e1.w * e2.dw) * S +
                        e1.w * e2.w * SDot;

                    if (e1.to == e2.to)
                    {
                        g_ADGraph->selfSoEdges[e1.to] += Real(2.0) * cross;
                        g_ADGraph->selfSoEdgesDot[e1.to] += Real(2.0) * crossDot;
                    }
                    else
                    {
                        g_ADGraph->soEdges[std::max(e1.to, e2.to)].Insert(
                            std::min(e1.to, e2.to),
                            cross);
                        g_ADGraph->soEdgesDot[std::max(e1.to, e2.to)].Insert(
                            std::min(e1.to, e2.to),
                            crossDot);
                    }
                }
            }

            // Create local second-order contribution and its directional derivative.
            const Real a = vertex.w;
            const Real aDot = vertex.wDot;

            if ((a != Real(0.0) || aDot != Real(0.0)) &&
                (vertex.soW != Real(0.0) || vertex.soWDot != Real(0.0)))
            {
                const Real create = a * vertex.soW;
                const Real createDot = aDot * vertex.soW + a * vertex.soWDot;

                if (e2.to == vid)
                {
                    g_ADGraph->selfSoEdges[e1.to] += create;
                    g_ADGraph->selfSoEdgesDot[e1.to] += createDot;
                }
                else if (e1.to == e2.to)
                {
                    g_ADGraph->selfSoEdges[e1.to] += Real(2.0) * create;
                    g_ADGraph->selfSoEdgesDot[e1.to] += Real(2.0) * createDot;
                }
                else
                {
                    g_ADGraph->soEdges[std::max(e1.to, e2.to)].Insert(
                        std::min(e1.to, e2.to),
                        create);
                    g_ADGraph->soEdgesDot[std::max(e1.to, e2.to)].Insert(
                        std::min(e1.to, e2.to),
                        createDot);
                }
            }

            // Propagate first-order adjoints and directional adjoints.
            if (a != Real(0.0) || aDot != Real(0.0))
            {
                vertex.w = Real(0.0);
                vertex.wDot = Real(0.0);

                g_ADGraph->vertices[e1.to].w += a * e1.w;
                g_ADGraph->vertices[e1.to].wDot +=
                    aDot * e1.w + a * e1.dw;

                if (e2.to != vid)
                {
                    g_ADGraph->vertices[e2.to].w += a * e2.w;
                    g_ADGraph->vertices[e2.to].wDot +=
                        aDot * e2.w + a * e2.dw;
                }
            }
        }
    }


    inline void PropagateAdjoint()
    {
        if (g_ADGraph->vertices.size() > g_ADGraph->soEdges.size())
        {
            g_ADGraph->soEdges.resize(g_ADGraph->vertices.size());
        }
        else
        {
            for (int i = 0; i < (int)g_ADGraph->soEdges.size(); i++)
            {
                g_ADGraph->soEdges[i].Clear();
            }
        }
        g_ADGraph->selfSoEdges.resize(g_ADGraph->vertices.size(), Real(0.0));
        // Any chance for SSE/AVX parallism?

        for (VertexId vid = g_ADGraph->vertices.size() - 1; vid > 0; vid--)
        {
            ADVertex &vertex = g_ADGraph->vertices[vid];
            ADEdge &e1 = vertex.e1;
            ADEdge &e2 = vertex.e2;
            if (e1.to == vid)
            {
                continue;
            }

            // Pushing
            BTree &btree = g_ADGraph->soEdges[vid];
            std::vector<BTNode>::iterator it;
            if (e2.to == vid)
            {
                for (it = btree.nodes.begin(); it != btree.nodes.end(); it++)
                {
                    ADEdge soEdge(it->key, it->val);
                    PushEdge(e1, soEdge);
                }
            }
            else
            {
                for (it = btree.nodes.begin(); it != btree.nodes.end(); it++)
                {
                    ADEdge soEdge(it->key, it->val);
                    PushEdge(e1, soEdge);
                    PushEdge(e2, soEdge);
                }
            }
            if (g_ADGraph->selfSoEdges[vid] != Real(0.0))
            {
                g_ADGraph->selfSoEdges[e1.to] += e1.w * e1.w * g_ADGraph->selfSoEdges[vid];
                if (e2.to != vid)
                {
                    g_ADGraph->selfSoEdges[e2.to] += e2.w * e2.w * g_ADGraph->selfSoEdges[vid];
                    if (e1.to == e2.to)
                    {
                        g_ADGraph->selfSoEdges[e2.to] += Real(2.0) * e1.w * e2.w * g_ADGraph->selfSoEdges[vid];
                    }
                    else
                    {
                        g_ADGraph->soEdges[std::max(e1.to, e2.to)].Insert(std::min(e1.to, e2.to),
                                                                          e1.w * e2.w * g_ADGraph->selfSoEdges[vid]);
                    }
                }
            }

            // release memory?

            Real a = vertex.w;
            if (a != Real(0.0))
            {
                // Creating
                if (vertex.soW != Real(0.0))
                {
                    if (e2.to == vid)
                    { // single-edge
                        g_ADGraph->selfSoEdges[e1.to] += a * vertex.soW;
                    }
                    else if (e1.to == e2.to)
                    {
                        g_ADGraph->selfSoEdges[e1.to] += Real(2.0) * a * vertex.soW;
                    }
                    else
                    {
                        g_ADGraph->soEdges[std::max(e1.to, e2.to)].Insert(std::min(e1.to, e2.to),
                                                                          a * vertex.soW);
                    }
                }
                // Adjoint
                vertex.w = Real(0.0);
                g_ADGraph->vertices[e1.to].w += a * e1.w;
                if (e2.to != vid)
                {
                    g_ADGraph->vertices[e2.to].w += a * e2.w;
                }
            }
        }
    }


    ////////////////// Quadra third-order extension API //////////////////

    typedef std::vector<std::vector<Real> > DenseMatrix;

    struct ValueGradientHessian
    {
        Real value = Real(0.0);
        std::vector<Real> gradient;
        DenseMatrix hessian;
    };

    template <typename Func>
    inline ValueGradientHessian evaluate_value_gradient_hessian(Func &&f,
                                                                const std::vector<Real> &x)
    {
        ADGraph graph;
        std::vector<AReal> ax;
        ax.reserve(x.size());
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            ax.push_back(AReal(x[i]));
        }

        AReal y = f(ax);
        SetAdjoint(y, Real(1.0));
        PropagateAdjoint();

        ValueGradientHessian out;
        out.value = y.val;
        out.gradient.resize(x.size());
        out.hessian.assign(x.size(), std::vector<Real>(x.size(), Real(0.0)));
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            out.gradient[i] = GetAdjoint(ax[i]);
            for (std::size_t j = 0; j < x.size(); ++j)
            {
                out.hessian[i][j] = GetAdjoint(ax[i], ax[j]);
            }
        }
        return out;
    }

    struct ThirdOrderScalar
    {
        Real val;
        Real d1;
        Real d2;
        Real d3;

        ThirdOrderScalar() : val(Real(0.0)), d1(Real(0.0)), d2(Real(0.0)), d3(Real(0.0)) {}
        ThirdOrderScalar(const Real v) : val(v), d1(Real(0.0)), d2(Real(0.0)), d3(Real(0.0)) {}
        ThirdOrderScalar(const Real v, const Real direction)
            : val(v), d1(direction), d2(Real(0.0)), d3(Real(0.0)) {}
    };

    struct DirectionalDerivatives3
    {
        Real value = Real(0.0);
        Real first = Real(0.0);   // df(x)[d]
        Real second = Real(0.0);  // d^T H(x) d
        Real third = Real(0.0);   // D^3 f(x)[d,d,d]
    };

    inline ThirdOrderScalar make_third_order_seed(const Real value, const Real direction)
    {
        return ThirdOrderScalar(value, direction);
    }

    inline ThirdOrderScalar unary_chain(const ThirdOrderScalar &x,
                                        const Real value,
                                        const Real fp,
                                        const Real fpp,
                                        const Real fppp)
    {
        ThirdOrderScalar y;
        y.val = value;
        y.d1 = fp * x.d1;
        y.d2 = fpp * x.d1 * x.d1 + fp * x.d2;
        y.d3 = fppp * x.d1 * x.d1 * x.d1 + Real(3.0) * fpp * x.d1 * x.d2 + fp * x.d3;
        return y;
    }

    inline ThirdOrderScalar operator+(const ThirdOrderScalar &a, const ThirdOrderScalar &b)
    {
        ThirdOrderScalar y;
        y.val = a.val + b.val;
        y.d1 = a.d1 + b.d1;
        y.d2 = a.d2 + b.d2;
        y.d3 = a.d3 + b.d3;
        return y;
    }
    inline ThirdOrderScalar operator+(const ThirdOrderScalar &a, const Real b) { return a + ThirdOrderScalar(b); }
    inline ThirdOrderScalar operator+(const Real a, const ThirdOrderScalar &b) { return ThirdOrderScalar(a) + b; }
    inline ThirdOrderScalar &operator+=(ThirdOrderScalar &a, const ThirdOrderScalar &b) { return (a = a + b); }
    inline ThirdOrderScalar &operator+=(ThirdOrderScalar &a, const Real b) { return (a = a + b); }

    inline ThirdOrderScalar operator-(const ThirdOrderScalar &a, const ThirdOrderScalar &b)
    {
        ThirdOrderScalar y;
        y.val = a.val - b.val;
        y.d1 = a.d1 - b.d1;
        y.d2 = a.d2 - b.d2;
        y.d3 = a.d3 - b.d3;
        return y;
    }
    inline ThirdOrderScalar operator-(const ThirdOrderScalar &a, const Real b) { return a - ThirdOrderScalar(b); }
    inline ThirdOrderScalar operator-(const Real a, const ThirdOrderScalar &b) { return ThirdOrderScalar(a) - b; }
    inline ThirdOrderScalar operator-(const ThirdOrderScalar &x)
    {
        ThirdOrderScalar y;
        y.val = -x.val;
        y.d1 = -x.d1;
        y.d2 = -x.d2;
        y.d3 = -x.d3;
        return y;
    }
    inline ThirdOrderScalar &operator-=(ThirdOrderScalar &a, const ThirdOrderScalar &b) { return (a = a - b); }
    inline ThirdOrderScalar &operator-=(ThirdOrderScalar &a, const Real b) { return (a = a - b); }

    inline ThirdOrderScalar operator*(const ThirdOrderScalar &a, const ThirdOrderScalar &b)
    {
        ThirdOrderScalar y;
        y.val = a.val * b.val;
        y.d1 = a.d1 * b.val + a.val * b.d1;
        y.d2 = a.d2 * b.val + Real(2.0) * a.d1 * b.d1 + a.val * b.d2;
        y.d3 = a.d3 * b.val + Real(3.0) * a.d2 * b.d1 + Real(3.0) * a.d1 * b.d2 + a.val * b.d3;
        return y;
    }
    inline ThirdOrderScalar operator*(const ThirdOrderScalar &a, const Real b)
    {
        ThirdOrderScalar y;
        y.val = a.val * b;
        y.d1 = a.d1 * b;
        y.d2 = a.d2 * b;
        y.d3 = a.d3 * b;
        return y;
    }
    inline ThirdOrderScalar operator*(const Real a, const ThirdOrderScalar &b) { return b * a; }
    inline ThirdOrderScalar &operator*=(ThirdOrderScalar &a, const ThirdOrderScalar &b) { return (a = a * b); }
    inline ThirdOrderScalar &operator*=(ThirdOrderScalar &a, const Real b) { return (a = a * b); }

    inline ThirdOrderScalar Inv(const ThirdOrderScalar &x)
    {
        const Real invX = Real(1.0) / x.val;
        const Real invX2 = invX * invX;
        const Real invX3 = invX2 * invX;
        const Real invX4 = invX3 * invX;
        return unary_chain(x, invX, -invX2, Real(2.0) * invX3, -Real(6.0) * invX4);
    }

    inline ThirdOrderScalar operator/(const ThirdOrderScalar &a, const ThirdOrderScalar &b) { return a * Inv(b); }
    inline ThirdOrderScalar operator/(const ThirdOrderScalar &a, const Real b) { return a * (Real(1.0) / b); }
    inline ThirdOrderScalar operator/(const Real a, const ThirdOrderScalar &b) { return a * Inv(b); }
    inline ThirdOrderScalar &operator/=(ThirdOrderScalar &a, const ThirdOrderScalar &b) { return (a = a / b); }
    inline ThirdOrderScalar &operator/=(ThirdOrderScalar &a, const Real b) { return (a = a / b); }

    inline ThirdOrderScalar square(const ThirdOrderScalar &x) { return x * x; }

    inline ThirdOrderScalar sqrt(const ThirdOrderScalar &x)
    {
        const Real sqrtX = std::sqrt(x.val);
        const Real invSqrtX = Real(1.0) / sqrtX;
        return unary_chain(x, sqrtX,
                           Real(0.5) * invSqrtX,
                           -Real(0.25) * invSqrtX / x.val,
                           Real(0.375) * invSqrtX / (x.val * x.val));
    }

    inline ThirdOrderScalar pow(const ThirdOrderScalar &x, const Real a)
    {
        return unary_chain(x,
                           std::pow(x.val, a),
                           a * std::pow(x.val, a - Real(1.0)),
                           a * (a - Real(1.0)) * std::pow(x.val, a - Real(2.0)),
                           a * (a - Real(1.0)) * (a - Real(2.0)) * std::pow(x.val, a - Real(3.0)));
    }

    inline ThirdOrderScalar exp(const ThirdOrderScalar &x)
    {
        const Real expX = std::exp(x.val);
        return unary_chain(x, expX, expX, expX, expX);
    }

    inline ThirdOrderScalar log(const ThirdOrderScalar &x)
    {
        const Real invX = Real(1.0) / x.val;
        return unary_chain(x, std::log(x.val), invX, -invX * invX, Real(2.0) * invX * invX * invX);
    }

    inline ThirdOrderScalar sin(const ThirdOrderScalar &x)
    {
        return unary_chain(x, std::sin(x.val), std::cos(x.val), -std::sin(x.val), -std::cos(x.val));
    }

    inline ThirdOrderScalar cos(const ThirdOrderScalar &x)
    {
        return unary_chain(x, std::cos(x.val), -std::sin(x.val), -std::cos(x.val), std::sin(x.val));
    }

    inline ThirdOrderScalar tan(const ThirdOrderScalar &x)
    {
        const Real tanX = std::tan(x.val);
        const Real secX = Real(1.0) / std::cos(x.val);
        const Real sec2X = secX * secX;
        return unary_chain(x, tanX,
                           sec2X,
                           Real(2.0) * tanX * sec2X,
                           Real(2.0) * sec2X * sec2X + Real(4.0) * tanX * tanX * sec2X);
    }

    inline ThirdOrderScalar asin(const ThirdOrderScalar &x)
    {
        const Real tmp = Real(1.0) / (Real(1.0) - x.val * x.val);
        const Real sqrtTmp = std::sqrt(tmp);
        return unary_chain(x, std::asin(x.val), sqrtTmp, x.val * sqrtTmp * tmp,
                           (Real(1.0) + Real(2.0) * x.val * x.val) * sqrtTmp * tmp * tmp);
    }

    inline ThirdOrderScalar acos(const ThirdOrderScalar &x)
    {
        const Real tmp = Real(1.0) / (Real(1.0) - x.val * x.val);
        const Real sqrtTmp = std::sqrt(tmp);
        return unary_chain(x, std::acos(x.val), -sqrtTmp, -x.val * sqrtTmp * tmp,
                           -(Real(1.0) + Real(2.0) * x.val * x.val) * sqrtTmp * tmp * tmp);
    }

    template <typename Func>
    inline DirectionalDerivatives3 evaluate_directional_derivatives3(Func &&f,
                                                                      const std::vector<Real> &x,
                                                                      const std::vector<Real> &direction)
    {
        if (x.size() != direction.size())
        {
            throw std::invalid_argument("evaluate_directional_derivatives3: x and direction must have the same length");
        }
        std::vector<ThirdOrderScalar> ax;
        ax.reserve(x.size());
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            ax.push_back(make_third_order_seed(x[i], direction[i]));
        }
        ThirdOrderScalar y = f(ax);
        DirectionalDerivatives3 out;
        out.value = y.val;
        out.first = y.d1;
        out.second = y.d2;
        out.third = y.d3;
        return out;
    }

    template <typename Func>
    inline Real third_directional_derivative(Func &&f,
                                             const std::vector<Real> &x,
                                             const std::vector<Real> &direction)
    {
        return evaluate_directional_derivatives3(std::forward<Func>(f), x, direction).third;
    }

    template <typename Func>
    inline DenseMatrix hessian_directional_derivative_central_difference(Func &&f,
                                                                         const std::vector<Real> &x,
                                                                         const std::vector<Real> &direction,
                                                                         Real eps = Real(1e-5))
    {
        if (x.size() != direction.size())
        {
            throw std::invalid_argument("hessian_directional_derivative_central_difference: x and direction must have the same length");
        }
        std::vector<Real> xp = x;
        std::vector<Real> xm = x;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            xp[i] += eps * direction[i];
            xm[i] -= eps * direction[i];
        }
        ValueGradientHessian hp = evaluate_value_gradient_hessian(f, xp);
        ValueGradientHessian hm = evaluate_value_gradient_hessian(f, xm);
        DenseMatrix out = hp.hessian;
        for (std::size_t i = 0; i < out.size(); ++i)
        {
            for (std::size_t j = 0; j < out[i].size(); ++j)
            {
                out[i][j] = (hp.hessian[i][j] - hm.hessian[i][j]) / (Real(2.0) * eps);
            }
        }
        return out;
    }


} // namespace had

#endif // HAD_QUADRA_H__


// ------------------------ base ------------------------
class CSR {
    constructor(n, edges) {
        this.start = Array(n + 1).fill(0);
        this.elist = Array(edges.length);

        // 各頂点の次数を計算
        edges.forEach(edge => {
            this.start[edge[0] + 1]++;
        });

        // 累積和を計算
        for (let i = 1; i <= n; i++) {
            this.start[i] += this.start[i - 1];
        }

        const counter = [...this.start];
        edges.forEach(edge => {
            this.elist[counter[edge[0]]++] = edge[1];
        });
    }

    N() {
        return this.start.length - 1;
    }

    deg(i) {
        return this.start[i + 1] - this.start[i];
    }

    begin(i) {
        return this.start[i];
    }

    end(i) {
        return this.start[i + 1];
    }

    at(i) {
        return this.elist[i];
    }

    static toCSR(graph) {
        const edges = [];
        graph.forEach((adjList, i) => {
            adjList.forEach(e => edges.push([i, e]));
        });
        return new CSR(graph.length, edges);
    }

    toVector() {
        const result = Array(this.N()).fill(null).map(() => []);
        for (let v = 0; v < this.N(); v++) {
            for (let i = this.begin(v); i < this.end(v); i++) {
                result[v].push(this.elist[i]);
            }
        }
        return result;
    }

    toEdgeList() {
        const result = [];
        for (let v = 0; v < this.N(); v++) {
            for (let i = this.begin(v); i < this.end(v); i++) {
                result.push([v, this.elist[i]]);
            }
        }
        return result;
    }
}

class edge {
    constructor(a, b, c) {
        this.s = a;
        this.t = b;
        this.w = c;
    }
};

class union_find {
    constructor(N) {
        this.sz = new Array(N);
        for (let i = 0; i < N; i++) this.sz[i] = -1;
    }

    find(i) {
        while (this.sz[i] >= 0) {
            i = this.sz[i];
        }
        return i;
    }

    size(i) {
        return -this.sz[this.find(i)];
    }

    same(a, b) {
        return this.find(a) == this.find(b);
    }

    unite(a, b) {
        a = this.find(a);
        b = this.find(b);
        if (a == b) return;
        if (this.sz[a] > this.sz[b]) {
            let t = a;
            a = b;
            b = t;
        }
        this.sz[a] += this.sz[b];
        this.sz[b] = a;
    }
};

class TECC {
    constructor(n) {
        this._n = n;
        this.edges = [];
    }

    numVertices() {
        return this._n;
    }

    addEdge(from, to) {
        this.edges.push([from, { to }]);
        this.edges.push([to, { to: from }]);
    }

    teccIDs() {
        const g = new CSR(this._n, this.edges);
        let nowOrd = 0;
        let groupNum = 0;
        const visited = [];
        const low = Array(this._n).fill(0);
        const ord = Array(this._n).fill(-1);
        const ids = Array(this._n).fill(0);

        const dfs = (v, p) => {
            low[v] = ord[v] = nowOrd++;
            visited.push(v);
            let revGuard = true;

            for (let i = g.start[v]; i < g.start[v + 1]; i++) {
                const to = g.elist[i].to;
                if (to === p && revGuard) {
                    revGuard = false;
                    continue;
                }
                if (ord[to] === -1) {
                    dfs(to, v);
                    low[v] = Math.min(low[v], low[to]);
                } else {
                    low[v] = Math.min(low[v], ord[to]);
                }
            }

            if (low[v] === ord[v]) {
                while (true) {
                    const u = visited.pop();
                    ord[u] = this._n;
                    ids[u] = groupNum;
                    if (u === v) break;
                }
                groupNum++;
            }
        };

        for (let i = 0; i < this._n; i++) {
            if (ord[i] === -1) dfs(i, -1);
        }

        ids.forEach((_, i) => {
            ids[i] = groupNum - 1 - ids[i];
        });

        return ids;
    }
}

class general_matching {
    constructor(n) {
        this.n = n;
        this.m = n + Math.floor(n / 2);
        this.mate = new Array(n).fill(-1);
        this.b = Array.from({ length: this.m }, () => []);
        this.p = new Array(this.m);
        this.d = new Array(this.m);
        this.bl = new Array(this.m);
        this.g = Array.from({ length: this.m }, () => new Array(this.m).fill(-1));
    }

    addEdge(u, v) {
        this.g[u][v] = u;
        this.g[v][u] = v;
    }

    match(u, v) {
        this.g[u][v] = this.g[v][u] = -1;
        this.mate[u] = v;
        this.mate[v] = u;
    }

    trace(x) {
        const vx = [];
        while (true) {
            while (this.bl[x] !== x) x = this.bl[x];
            if (vx.length && vx[vx.length - 1] === x) break;
            vx.push(x);
            x = this.p[x];
        }
        return vx;
    }

    contract(c, x, y, vx, vy) {
        this.b[c] = [];
        let r = vx[vx.length - 1];
        while (vx.length && vy.length && vx[vx.length - 1] === vy[vy.length - 1]) {
            r = vx.pop();
            vy.pop();
        }
        this.b[c].push(r);
        this.b[c].push(...vx.reverse());
        this.b[c].push(...vy);

        for (let i = 0; i <= c; i++) {
            this.g[c][i] = this.g[i][c] = -1;
        }

        for (let z of this.b[c]) {
            this.bl[z] = c;
            for (let i = 0; i < c; i++) {
                if (this.g[z][i] !== -1) {
                    this.g[c][i] = z;
                    this.g[i][c] = this.g[i][z];
                }
            }
        }
    }

    lift(vx) {
        const A = [];
        while (vx.length >= 2) {
            let z = vx.pop();
            if (z < this.n) {
                A.push(z);
                continue;
            }
            let w = vx[vx.length - 1];
            let i = (A.length % 2 === 0 ? this.b[z].indexOf(this.g[z][w]) : 0);
            let j = (A.length % 2 === 1 ? this.b[z].indexOf(this.g[z][A[A.length - 1]]) : 0);
            let k = this.b[z].length;
            let dif = (A.length % 2 === 0 ? i % 2 === 1 : j % 2 === 0) ? 1 : k - 1;

            while (i !== j) {
                vx.push(this.b[z][i]);
                i = (i + dif) % k;
            }
            vx.push(this.b[z][i]);
        }
        return A;
    }

    solve() {
        for (let ans = 0; ; ans++) {
            this.d.fill(0);
            const Q = [];
            for (let i = 0; i < this.m; i++) this.bl[i] = i;
            for (let i = 0; i < this.n; i++) {
                if (this.mate[i] === -1) {
                    Q.push(i);
                    this.p[i] = i;
                    this.d[i] = 1;
                }
            }
            let c = this.n;
            let aug = false;

            while (Q.length && !aug) {
                let x = Q.shift();
                if (this.bl[x] !== x) continue;

                for (let y = 0; y < c; y++) {
                    if (this.bl[y] === y && this.g[x][y] !== -1) {
                        if (this.d[y] === 0) {
                            this.p[y] = x;
                            this.d[y] = 2;
                            this.p[this.mate[y]] = y;
                            this.d[this.mate[y]] = 1;
                            Q.push(this.mate[y]);
                        } else if (this.d[y] === 1) {
                            const vx = this.trace(x);
                            const vy = this.trace(y);

                            if (vx[vx.length - 1] === vy[vy.length - 1]) {
                                this.contract(c, x, y, vx, vy);
                                Q.push(c);
                                this.p[c] = this.p[this.b[c][0]];
                                this.d[c] = 1;
                                c++;
                            } else {
                                aug = true;
                                vx.unshift(y);
                                vy.unshift(x);
                                const A = this.lift(vx);
                                const B = this.lift(vy);
                                A.push(...B.reverse());

                                for (let i = 0; i < A.length; i += 2) {
                                    this.match(A[i], A[i + 1]);
                                    if (i + 2 < A.length) this.addEdge(A[i + 1], A[i + 2]);
                                }
                            }
                            break;
                        }
                    }
                }
            }
            if (!aug) return ans;
        }
    }

    // マッチ先(存在しない場合-1)
    getMate(i) {
        return this.mate[i];
    }
}

class SCC {
    constructor(n = 0) {
        this._n = n;
        this.edges = [];
    }

    numVertices() {
        return this._n;
    }

    addEdge(from, to) {
        this.edges.push([from, { to }]);
    }

    sccIds() {
        const g = new CSR(this._n, this.edges);
        let nowOrd = 0;
        let groupNum = 0;
        const visited = [];
        const low = new Array(this._n).fill(0);
        const ord = new Array(this._n).fill(-1);
        const ids = new Array(this._n);

        const dfs = (v) => {
            low[v] = ord[v] = nowOrd++;
            visited.push(v);

            for (let i = g.start[v]; i < g.start[v + 1]; i++) {
                const to = g.elist[i].to;
                if (ord[to] === -1) {
                    dfs(to);
                    low[v] = Math.min(low[v], low[to]);
                } else {
                    low[v] = Math.min(low[v], ord[to]);
                }
            }

            if (low[v] === ord[v]) {
                while (true) {
                    const u = visited.pop();
                    ord[u] = this._n;
                    ids[u] = groupNum;
                    if (u === v) break;
                }
                groupNum++;
            }
        };

        for (let i = 0; i < this._n; i++) {
            if (ord[i] === -1) dfs(i);
        }

        for (let i = 0; i < ids.length; i++) {
            ids[i] = groupNum - 1 - ids[i];
        }

        return ids;
    }
}

class Dijkstra {
    constructor(N) {
        this.graph = new Array();
        for (let i = 0; i < N; i++) {
            this.graph.push(new Array());
        }
        this.N = N;
        this.dist = [];
        this.par = [];
        this._s = -1;
        this.INF = Number.MAX_SAFE_INTEGER / 2;
    }

    addEdge(s, t, w) {
        this.graph[s].push(new edge(s, t, w));
    }

    build(start) {
        this._s = start;
        const N = this.N;
        const INF = this.INF;
        // 距離と親ノードの初期化
        if (this.dist.length === 0) {
            this.dist = Array(N).fill(INF);
            this.par = Array(N).fill(-1);
        } else {
            this.dist.fill(INF);
            this.par.fill(-1);
        }
        const priorityQueue = [];
        const enqueue = (weight, vertex) => {
            priorityQueue.push({ weight, vertex });
            priorityQueue.sort((a, b) => a.weight - b.weight);
        };
        const dequeue = () => priorityQueue.shift();
        this.dist[start] = 0;
        enqueue(0, start);
        while (priorityQueue.length > 0) {
            const { weight: w, vertex: v } = dequeue();
            if (this.dist[v] < w) continue;
            for (let i = 0; i < this.graph[v].length; i++) {
                const edge = this.graph[v][i];
                const to = edge.t;
                const d = this.dist[v] + edge.w;
                if (this.dist[to] > d) {
                    this.dist[to] = d;
                    this.par[to] = v;
                    enqueue(d, to);
                }
            }
        }
    }
}

class Dinic {
    constructor(n = 0) {
        this.N = n; // 頂点数
        this.pos = []; // 辺の情報を保持
        this.G = Array.from({ length: n }, () => []); // 隣接リスト
    }

    // 頂点を追加
    addVertex() {
        this.G.push([]);
        return this.N++;
    }

    // 辺を追加
    addEdge(s, t, cap) {
        const m = this.pos.length;
        this.pos.push([s, this.G[s].length]);
        this.G[s].push({ t, rev: this.G[t].length, cap });
        this.G[t].push({ t: s, rev: this.G[s].length - 1, cap: 0 });
        return m;
    }

    // 辺情報を取得
    getEdge(i) {
        const [s, idx] = this.pos[i];
        const edge = this.G[s][idx];
        const revEdge = this.G[edge.t][edge.rev];
        return {
            s: s,
            t: edge.t,
            cap: edge.cap + revEdge.cap,
            flow: revEdge.cap,
        };
    }

    // 全ての辺情報を取得
    edges() {
        return this.pos.map((_, i) => this.getEdge(i));
    }

    // 辺を変更
    changeEdge(i, newCap, newFlow) {
        const [s, idx] = this.pos[i];
        const edge = this.G[s][idx];
        const revEdge = this.G[edge.t][edge.rev];
        edge.cap = newCap - newFlow;
        revEdge.cap = newFlow;
    }

    // 最大流を計算
    maxFlow(s, t, flowLimit = Infinity) {
        const level = Array(this.N).fill(-1);
        const iter = Array(this.N).fill(0);

        const bfs = () => {
            level.fill(-1);
            level[s] = 0;
            const que = [s];
            let l = 0;
            while (l < que.length) {
                const v = que[l++];
                for (const e of this.G[v]) {
                    if (e.cap > 0 && level[e.t] === -1) {
                        level[e.t] = level[v] + 1;
                        if (e.t === t) return;
                        que.push(e.t);
                    }
                }
            }
        };

        const dfs = (v, up) => {
            if (v === s) return up;
            let res = 0;
            const levelV = level[v];
            for (; iter[v] < this.G[v].length; iter[v]++) {
                const e = this.G[v][iter[v]];
                if (levelV <= level[e.t] || this.G[e.t][e.rev].cap === 0) continue;
                const d = dfs(e.t, Math.min(up - res, this.G[e.t][e.rev].cap));
                if (d <= 0) continue;
                e.cap += d;
                this.G[e.t][e.rev].cap -= d;
                res += d;
                if (res === up) break;
            }
            return res;
        };

        let flow = 0;
        while (flow < flowLimit) {
            bfs();
            if (level[t] === -1) break;
            iter.fill(0);
            while (flow < flowLimit) {
                const f = dfs(t, flowLimit - flow);
                if (f === 0) break;
                flow += f;
            }
        }
        return flow;
    }

    // s側の頂点かどうか (maxFlowを先に呼び出す必要あり)
    minCut(s) {
        const visited = Array(this.N).fill(false);
        const que = [s];
        let l = 0;
        while (l < que.length) {
            const p = que[l++];
            visited[p] = true;
            for (const e of this.G[p]) {
                if (e.cap > 0 && !visited[e.t]) {
                    visited[e.t] = true;
                    que.push(e.t);
                }
            }
        }
        return visited;
    }
}
// ------------------------ random graph ------------------------

// O(Mlog(M))
// O(Mlog^2(M)) (connected)
function random_undirected_graph(N, M, use_selfedge, use_multiedge, connected, minw = 1, maxw = 1) {
    if (N == 0 || (N == 1 && !use_selfedge)) return [];
    if (connected && M < N - 1) return [];
    let uf = new union_find(N);
    let res = [];
    if (!use_multiedge) {
        let maxM = (N * (N - 1)) / 2 + (use_selfedge ? N : 0);
        if (maxM < M) return [];
        let used = new Set();
        for (let i = 0; i < M; i++) {
            let a, b;
            while (true) {
                a = random_int(N - 1);
                b = random_int(N - 1);
                if (a > b) {
                    let tmp = a;
                    a = b;
                    b = tmp;
                }
                if (!use_selfedge && a == b) continue;
                if (connected && uf.size(0) != N && uf.same(a, b)) continue;
                let c = a * N + b;
                if (used.has(c)) continue;
                used.add(c);
                res.push(new edge(a, b, random_int_range(minw, maxw)));
                uf.unite(a, b);
                break;
            }
        }
    } else {
        for (let i = 0; i < M; i++) {
            let a, b;
            while (true) {
                a = random_int(N - 1);
                b = random_int(N - 1);
                if (a > b) {
                    let tmp = a;
                    a = b;
                    b = tmp;
                }
                if (!use_selfedge && a == b) continue;
                if (connected && uf.size(0) != N && uf.same(a, b)) continue;
                res.push(new edge(a, b, random_int_range(minw, maxw)));
                uf.unite(a, b);
                break;
            }
        }
    }
    return res;
}

// O(MlogM)
function random_directed_graph(N, M, use_selfedge, use_multiedge, minw = 1, maxw = 1) {
    if (N == 0 || (N == 1 && !use_selfedge)) return [];
    let res = [];
    if (!use_multiedge) {
        let maxM = N * (N - 1) + (use_selfedge ? N : 0);
        if (maxM < M) return [];
        let used = new Set();
        for (let i = 0; i < M; i++) {
            let a, b;
            while (true) {
                a = random_int(N - 1);
                b = random_int(N - 1);
                if (!use_selfedge && a == b) continue;
                let c = a * N + b;
                if (used.has(c)) continue;
                used.add(c);
                res.push(new edge(a, b, random_int_range(minw, maxw)));
                break;
            }
        }
    } else {
        for (let i = 0; i < M; i++) {
            let a, b;
            while (true) {
                a = random_int(N - 1);
                b = random_int(N - 1);
                if (!use_selfedge && a == b) continue;
                res.push(new edge(a, b, random_int_range(minw, maxw)));
                break;
            }
        }
    }
    return res;
}


// ------------------------ calc graph info ------------------------

// O(N^2)
function calc_complement_graph(N, E) {
    let used = [];
    for (let i = 0; i < N; i++) {
        used.push(new Array(N).fill(0));
    }
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        used[e.s][e.t] = 1;
        used[e.t][e.s] = 1;
    }
    let _E = []
    for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
            if (!used[i][j]) {
                _E.push(new edge(i, j, 1, -1));
            }
        }
    }
    return _E;
}

// O(NlogN)
function calc_connected_component(N, E) {
    let uf = new union_find(N);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        uf.unite(e.s, e.t);
    }
    let res = new Array(N);
    for (let i = 0; i < N; i++) {
        let r = uf.find(i);
        res[i] = r;
    }
    return res;
}

// O(N + M)
function calc_two_edge_connected_component(N, E) {
    let g = new TECC(N);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        g.addEdge(e.s, e.t);
    }
    return g.teccIDs();
}

// O((N + M)log(N))
function calc_bipartite_graph(N, E) {
    let uf = new union_find(2 * N);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        uf.unite(e.s, e.t + N);
        uf.unite(e.s + N, e.t);
    }
    for (let i = 0; i < N; i++) {
        if (uf.same(i, i + N)) {
            return [false, []];
        }
    }
    let cc = new Array(2 * N, 0);
    for (let i = 0; i < 2 * N; i++) {
        let r = uf.find(i);
        cc[r].push(i);
    }
    let ans = new Array(N).fill(-1);
    for (let i = 0; i < 2 * N; i++) {
        if (cc[i].length == 0) continue;
        let t = (cc[i][0] >= 0 ? cc[i][0] - N : cc[i][0]);
        if (ans[t] != -1) continue;
        for (let j = 0; j < cc[i].length; j++) {
            let v = cc[i][j];
            if (v >= N) ans[v - N] = 1;
            else ans[v] = 0;
        }
    }
    return [true, ans];
}

// O(2^N * N)
function calc_chromatic_number(N, E) {
    const n = N;
    const m = 1 << n;
    if (n === 0) return [];
    const z = Array(m).fill(0); // 1色で塗れるか
    const g = Array(n).fill(0);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        if (e.s == e.t) return [];
        g[e.s] |= 1 << e.t;
        g[e.t] |= 1 << e.s;
    }
    for (let i = 0; i < m; i++) {
        let ok = true;
        for (let j = 0; j < n; j++) {
            if (((i >> j) & 1) && (i & g[j])) {
                ok = false;
                break;
            }
        }
        z[i] = ok ? 1 : 0;
    }
    const zeta = (A) => {
        const k = A.length;
        for (let i = 1; i < k; i <<= 1) {
            for (let j = 0; j < k; j++) {
                if (!(j & i)) {
                    A[j | i] += A[j];
                }
            }
        }
    };
    const mobius = (A) => {
        const k = A.length;
        for (let i = 1; i < k; i <<= 1) {
            for (let j = 0; j < k; j++) {
                if (!(j & i)) {
                    A[j | i] -= A[j];
                }
            }
        }
    };
    const dp = [];
    dp.push([...z]);
    for (let i = 0;; i++) {
        if (!dp[dp.length - 1][(1 << (n - i)) - 1]) {
            const k = 1 << (n - i - 1);
            const x = Array(k).fill(0);
            const y = Array(k).fill(0);
            for (let j = 0; j < (1 << n); j++) {
                x[j >> (i + 1)] |= z[j];
            }
            for (let j = 1; j < 2 * k; j += 2) {
                y[j >> 1] |= dp[dp.length - 1][j] ? 1 : 0;
            }
            zeta(x);
            zeta(y);
            for (let j = 0; j < k; j++) {
                x[j] *= y[j];
            }
            mobius(x);
            dp.push(x);
        } else {
            const res = Array(i + 1).fill().map(() => []);
            let unused = (1 << n) - 1;
            for (; i >= 0; i--) {
                const mask = (1 << i) - 1;
                for (let j = 1; j < (1 << n); j++) {
                    if (!z[j] || (j & unused) !== j) continue;
                    const f = unused ^ j;
                    if (i) {
                        if ((f & mask) !== mask || !dp[i - 1][f >> (i - 1)]) continue;
                    } else {
                        if (f) continue;
                    }
                    for (let k = 0; k < n; k++) {
                        if ((j >> k) & 1) {
                            res[i].push(k);
                        }
                    }
                    unused = f;
                    break;
                }
            }
            return res;
        }
    }
}

// O(N^3)
function calc_general_matching(N, E) {
    let g = new general_matching(N);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        g.addEdge(e.s, e.t);
    }
    let m = g.solve();
    return [m, g.mate];
}

// O(N + M)
function calc_strongly_connected_component(N, E) {
    let g = new SCC(N);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        g.addEdge(e.s, e.t);
    }
    return g.sccIds();
}

// O((N + M)log(N))
function calc_dijkstra(N, E, s) {
    let g = new Dijkstra(N);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        g.addEdge(e.s, e.t, e.w);
    }
    g.build(s);
    return [g.dist, g.par];
}

// O(N^2M)
function calc_max_flow_min_cut(N, E, s, t) {
    let g = new Dinic(N);
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        g.addEdge(e.s, e.t, e.w);
    }
    let f = g.maxFlow(s, t);
    let _E = [];
    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        let w = g.getEdge(i).flow;
        if (w != 0) _E.push(new edge(e.s, e.t, w));
    }
    return [f, _E, g.minCut(s)];
}
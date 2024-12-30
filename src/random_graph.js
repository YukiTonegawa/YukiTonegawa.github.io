import * as rnd from './random.js'

export class graph {
    constructor(N) {
        this.N = N;
        this.E = []
        this.is_directed = false;
    }
    add_edge(a, b, c) {
        this.E.push([a, b, c])
    }
};

export class union_find {
    constructor(N) {
        this.sz = new Array(N);
        for (let i = 0; i < N; i++) this.sz[i] = -1;
    }

    find(i) {
        if (this.sz[i] < 0) return i;
        return this.sz[i] = this.find(this.sz[i]);
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
        if (this.sz[a] > this.sz[b]) {
            let t = a;
            a = b;
            b = t;
        }
        this.sz[a] += this.sz[b];
        this.sz[b] = a;
    }
};

// 単純無向グラフ
// 0 <= N, 0 <= M <= NC2
// O(N^2)
export function random_undirected_simple_graph(N, M) {
    let maxM = N * (N - 1) / 2;
    let res = new graph(N);
    let flag = new Array(maxM);
    for (let i = 0; i < maxM; i++) flag[i] = (i < M);
    flag = rnd.random_shuffle(flag);
    let pos = 0;
    for (let a = 0; a < N; a++) {
        for (let b = a + 1; b < N; b++) {
            if (flag[pos++]) {
                res.add_edge(a, b, 1);
            }
        }
    }
    return res;
}

// O(N^2α(N))
export function random_spanning_tree(N) {
    let maxM = N * (N - 1) / 2;
    let res = new graph(N);
    let left = [];
    for (let a = 0; a < N; a++) {
        for (let b = a + 1; b < N; b++) {
            left.push(a);
        }
    }
    let uf = new union_find(N);
    let order = rnd.random_permutation(maxM);
    for (let p in order) {
        let a = left[p];
        let b = p % (N - a) + a;
        if (!uf.same(a, b)) {
            res.add_edge(a, b, 1);
            uf.unite(a, b);
            if (uf.size(0) == N) break;
        }
    }
    return res;
}

// 単純無向連結グラフ M >= N - 1
// 0 <= N, N-1 <= M <= NC2
// O(N^2α(N))
export function random_undirected_simple_connected_graph(N, M) {
    let minM = N - 1, maxM = N * (N - 1) / 2;
    let res = new graph(N);
    let left = [];
    for (let a = 0; a < N; a++) {
        for (let b = a + 1; b < N; b++) {
            left.push(a);
        }
    }
    let uf = new union_find(N);
    let order = rnd.random_permutation(maxM);
    for (let p in order) {
        let a = left[p - 1];
        let b = (p - 1) % (N - a) + a;
        if (!uf.same(a, b)) {
            res.add_edge(a, b, 1);
            uf.unite(a, b);
            left[p - 1] = -1; 
            if (uf.size(0) == N) break;
        }
    }
    order = rnd.random_shuffle(order);
    let ecnt = minM;
    for (let p in order) {
        let a = left[p - 1];
        if (a == -1) continue;
        let b = (p - 1) % (N - a) + a;
        res.add_edge(a, b, 1);
        ecnt++;
        if (ecnt == M) return res;
    }
    return res;
}

// pが大きいと直線に近く、小さいと完全二分木、ウニに近い
export function random_tree(N, p) {
    //
}
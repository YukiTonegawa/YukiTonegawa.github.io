// V[l, r)
export function get_subarray(V, l, r) {
    return V.slice(l, r);
}

export function calc_sum(V) {
    let res = 0;
    for (let x in V) {
        res += x;
    }
    return res;
}

export function calc_min(V) {
    return Math.min(...V);
}

export function calc_max(V) {
    return Math.max(...V);
}

export function calc_avg(V) {
    return calc_sum(V) / V.length;
}

export function gcd(a, b) {
    if (a < b) return gcd(b, a);
    let tmp;
    while (b) {
        tmp = a % b;
        a = b;
        b = tmp;
    }
    return a;
}

export function calc_gcd(V) {
    let res = 0;
    for (let x in V) {
        res = gcd(res, x);
    }
    return res;
}

export function calc_integral(V) {
    let N = V.length;
    for (let i = 1; i < N; i++) {
        V[i] += V[i - 1];
    }
    return V;
}

export function calc_diff(V) {
    let N = V.length;
    for (let i = N - 1; i >= 0; i--) {
        V[i] -= V[i - 1];
    }
    return V;
}

// O(N^2)
// indices of lis
export function calc_lis(V) {
    let N = V.length;
    if (N == 0) return [];
    let len = [];
    let prev = [];
    for (let i = 0; i < N; i++) {
        let l = 1;
        let p = -1;
        for (let j = 0; j < i; j++) {
            if (V[j] < V[i] && l < len[j] + 1) {
                l = len[j] + 1;
                p = j;
            }
        }
        len.push(l);
        prev.push(p);
    }
    let Mlen = Math.max(...len);
    let res = [];
    for (let i = 0; i < N; i++) {
        if (len[i] == Mlen) {
            let res = [];
            while (i >= 0) {
                res.push(i);
                i = prev[i];
            }
            return res.reverse();
        }
    }
    return [];
}

// O(N^2)
export function calc_inversion(V) {
    let N = V.length;
    let res = 0;
    for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
            res += (V[i] > V[j]);
        }
    }
    return res;
}

export function calc_mex(V) {
    let N = V.length;
    // 0 <= mex <= Nなので全探索
    let contained = [];
    for (let i = 0; i <= N; i++) contained.push(0);
    for (let x in V) {
        if (0 <= x && x <= N) {
            contained[x] = 1;
        }
    }
    for (let i = 0; i <= N; i++) {
        if (!contained[i]) {
            return i;
        }
    }
    return N;
}

// 和が最大の連続部分列[l, r)(空の列も和0として許容する) 　{l, r, sum}を返す
export function calc_maximum_subarray(V) {
    let N = V.length;
    let Msum = 0;
    let Ml = 0, Mr = 0;
    for (let l = 0; l < N; l++) {
        let r = l;
        let s = 0;
        while (r < N) {
            s += V[r];
            r++;
            if (Msum < s) {
                Msum = s;
                Ml = l, Mr = r;
            }
        }
    }
    return [l, r, Msum];
}

// 連続とは限らない部分列の個数(空の列も数える)
export function calc_number_of_subsequence(V) {
    let N = V.length;
    let dp = [1];
    let mp = new Map();
    for (let i = 0; i < N; i++) {
        let x = dp[i] * 2;
        if (mp.has(V[i])) x -= dp[mp.get(V[i])];
        dp.push(x);
        mp.set(V[i], i);
    }
    return dp[N];
}

export class cartesian_tree {
    find_same(i) {
        i = this.same[i];
        if (this.p[i] == -1 || this.val[p[i]] != this.val[i]) return i;
        return this.same[i] = this.find_same(this.p[i]);
    }

    constructor(V) {
        let N = V.length;
        this.val = V;
        this.l = Array(N);
        this.r = Array(N);
        this.p = Array(N);
        this.same = Array(N);
        for (let i = 0; i < N; i++) {
            this.r[i] = N;
            this.p[i] = -1;
            this.same[i] = i;
        }
        let st = [];
        for (let i = 0; i < N; i++) {
            let last = -1;
            while (st.length > 0 && V[i] < V[st[st.length - 1]]) {
                last = st.pop();
                this.r[last] = i;
            }
            if (last != -1) this.p[last] = i;
            const len = st.length;
            if (len > 0) {
                this.p[i] = st[len - 1];
                this.l[i] = st[len - 1] + 1;
            } else {
                this.l[i] = 0;
            }
            st.push(i);
        }
    }

    // 親ノード, 根の場合は-1を返す
    parent(i) { 
        return this.p[i];
    }
    
    // 部分木サイズ
    size_subtree(i) {
        return this.r[i] - this.l[i];
    }
  
    // v[i]の部分木を表す半開区間[l, r)
    segment_subtree(i) {
        return this.l[i], this.r[i];
    }
  
    // v[i]が最小であるような最大の半開区間[l, r)
    // 値がユニークなら segment_subtree = segment_minimum
    segment_minumum(i) {
        return segment_subtree(this.find_same(i));
    }
};

// ヒストグラム中の最大長方形の[l, r), 面積を返す
export function calc_largest_rectangle_histogram(V) {
    let N = V.length;
    let T = new cartesian_tree(V);
    let Ml = 0, Mr = 0, Ms = 0;
    for (let i = 0; i < N; i++) {
        l, r = T.segment_subtree(i);
        if ((r - l) * V[i] > Ms) {
            Ml = l;
            Mr = r;
            Ms = (r - l) * V[i];
        }
    }
    return Ml, Mr, Ms;
}

// precursive
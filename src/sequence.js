function calc_sum(V) {
    let res = 0;
    for (let i = 0; i < V.length; i++) {
        res += V[i];
    }
    return res;
}

function calc_min(V) {
    return Math.min(...V);
}

function calc_max(V) {
    return Math.max(...V);
}

function calc_avg(V) {
    return calc_sum(V) / V.length;
}

function gcd(a, b) {
    if (a < b) return gcd(b, a);
    let tmp;
    while (b) {
        tmp = a % b;
        a = b;
        b = tmp;
    }
    return a;
}

function calc_gcd(V) {
    let res = 0;
    for (let i = 0; i < V.length; i++) {
        res = gcd(res, Math.abs(V[i]));
    }
    return res;
}

function calc_cumsum(V, op) {
    let N = V.length;
    let res = new Array(N);
    for (let i = 0; i < N; i++) {
        res[i] = V[i];
        if (i) res[i] = op(res[i - 1], res[i]);
    }
    return res;
}

function calc_diff(V) {
    let N = V.length;
    let res = new Array(N);
    for (let i = N - 1; i >= 0; i--) {
        res[i] = V[i] - (i == 0 ? 0 : V[i - 1]);
    }
    return res;
}

// O(N^2)
// indices of lis
function calc_lis(V) {
    let N = V.length;
    if (N == 0) return [];
    let len = [];
    let prev = [];
    for (let i = 0; i < N; i++) {
        let l = 1, p = -1;
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

function calc_lis_val(V) {
    let idx = calc_lis(V);
    let res = [];
    for (let i = 0; i < idx.length; i++) {
        res.push(V[idx[i]]);
    }
    return res;
}

// O(N^2)
function calc_inversion(V) {
    let N = V.length;
    let res = 0;
    for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
            res += (V[i] > V[j]);
        }
    }
    return res;
}

function calc_mex(V) {
    let N = V.length;
    // 0 <= mex <= Nなので全探索
    let contained = [];
    for (let i = 0; i <= N; i++) contained.push(0);
    for (let i = 0; i < N; i++) {
        if (0 <= V[i] && V[i] <= N) {
            contained[V[i]] = 1;
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
function calc_maximum_subarray(V) {
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
    return [Ml, Mr, Msum]
}

function calc_maximum_subarray_val(V) {
    let [l, r, s] = calc_maximum_subarray(V);
    let res = [];
    for (let i = l; i < r; i++) {
        res.push(V[i]);
    }
    return res;
}

// 連続とは限らない部分列の個数(空の列も数える)
function calc_number_of_subsequence(V) {
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

class cartesian_tree {
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
        return [this.l[i], this.r[i]];
    }
  
    // v[i]が最小であるような最大の半開区間[l, r)
    // 値がユニークなら segment_subtree = segment_minimum
    segment_minumum(i) {
        return segment_subtree(this.find_same(i));
    }
};

// ヒストグラム中の最大長方形の[l, r), 高さ, 面積を返す
function calc_largest_rectangle_histogram(V) {
    let N = V.length;
    let T = new cartesian_tree(V);
    let Ml = 0, Mr = 0, Mh = 0, Ms = 0;
    for (let i = 0; i < N; i++) {
        let [l, r] = T.segment_subtree(i);
        if ((r - l) * V[i] > Ms) {
            Ml = l;
            Mr = r;
            Mh = V[i];
            Ms = (r - l) * V[i];
        }
    }
    return [Ml, Mr, Mh, Ms];
}


function safeMod(x, m) {
    x = BigInt(x);
    m = BigInt(m);
    x %= m;
    if (x < 0n) x += m;
    return x;
}

function gcd(a, b) {
    a = BigInt(Math.abs(Number(a)));
    b = BigInt(Math.abs(Number(b)));
    while (b) {
        [a, b] = [b, a % b];
    }
    return a;
}

function invGcd(a, b) {
    a = safeMod(a, b);
    if (a === 0n) return { gcd: b, inv: 0n };
    let s = b, t = a;
    let m0 = 0n, m1 = 1n;
    while (t) {
        const u = s / t;
        [s, t] = [t, s - t * u];
        [m0, m1] = [m1, m0 - m1 * u];
    }
    if (m0 < 0n) m0 += b / s;
    return { gcd: s, inv: m0 };
}

function addMod(a, b, mod) {
    a = BigInt(a);
    b = BigInt(b);
    mod = BigInt(mod);
    return (a + b >= mod ? a + b - mod : a + b);
}

function subMod(a, b, mod) {
    a = BigInt(a);
    b = BigInt(b);
    mod = BigInt(mod);
    return (a < b ? a + (mod - b) : a - b);
}

function mulMod(a, b, mod) {
    a = BigInt(a);
    b = BigInt(b);
    mod = BigInt(mod);
    return (a * b) % mod;
}

function divMod(a, b, mod) {
    const { gcd, inv } = invGcd(b, mod);
    if (gcd !== 1n) throw new Error("Division not possible");
    return mulMod(a, inv, mod);
}

function gaussianElimination(v, mod) {
    const n = v.length, m = (n === 0 ? 0 : v[0].length);
    let row = 0;
    for (let i = 0; i < m && row < n; i++) {
        let r = -1;
        for (let j = row; j < n; j++) {
            if (v[j][i] !== 0n) {
                r = j;
                break;
            }
        }
        if (r === -1) continue;
        if (r !== row) [v[r], v[row]] = [v[row], v[r]];
        for (let j = row + 1; j < n; j++) {
            if (v[j][i] === 0n) continue;
            const x = divMod(v[j][i], v[row][i], mod);
            for (let k = i; k < m; k++) {
                v[j][k] = subMod(v[j][k], mulMod(x, v[row][k], mod), mod);
            }
        }
        row++;
    }
    return v;
}

function systemOfLinearEquations(A, b, mod) {
    const n = A.length, m = (n === 0 ? 0 : A[0].length) + 1;
    for (let i = 0; i < n; i++) {
        A[i].push(BigInt(b[i]));
    }
    A = gaussianElimination(A, mod);
    const rank = (() => {
        let cnt = 0;
        for (let i = 0; i < n; i++, cnt++) {
            if (!A[i].some(v => v !== 0n)) break;
        }
        return cnt;
    })();

    const col = Array(rank).fill(-1);
    for (let i = 0; i < rank; i++) {
        let inv = 0n;
        let f = false;
        for (let j = i; j < m; j++) {
            if (A[i][j] === 0n) continue;
            if (j === m - 1 && !f) {
                return { vec: [], sol: [] };
            }
            if (!f) {
                inv = invGcd(A[i][j], mod).inv;
                col[i] = j;
                f = true;
            }
            A[i][j] = mulMod(A[i][j], inv, mod);
        }
    }

    const d = m - 1 - rank, v = m - 1;
    const sol = Array(v).fill(0n);
    for (let i = rank - 1; i >= 0; i--) {
        const idx = col[i];
        sol[idx] = A[i][v];
        for (let j = idx + 1; j < v; j++) {
            sol[idx] = subMod(sol[idx], mulMod(sol[j], A[i][j], mod), mod);
        }
    }

    const vec = Array.from({ length: d }, () => Array(v).fill(0n));
    const notCol = Array(v).fill(true);
    for (const c of col) notCol[c] = false;
    for (let i = 0, j = 0; i < v; i++) {
        if (notCol[i]) vec[j++][i] = 1n;
    }

    for (let i = rank - 1; i >= 0; i--) {
        const c = col[i];
        for (let k = 0; k < d; k++) {
            for (let j = c + 1; j < v; j++) {
                vec[k][c] = subMod(vec[k][c], mulMod(vec[k][j], A[i][j], mod), mod);
            }
        }
    }

    return { vec, sol };
}

function findPRecursive(_A, _mod, r, d) {
    let mod = BigInt(_mod);
    if (r < 0 || d < 0 || mod <= 0) throw new Error("Invalid parameters");
    let A = [];
    const N = _A.length;
    for (let i = 0; i < N; i++) {
        let x = BigInt(_A[i]) % mod;
        if (x < 0) x += mod;
        A.push(x);
    }
    const m = (r + 1) * (d + 1);
    if (N < m + r) return [];
    const M = Array.from({ length: N - r }, () => Array(m).fill(0n));
    for (let k = 0; k < N - r; k++) {
        for (let i = 0; i <= r; i++) {
            let x = BigInt(A[k + i]);
            for (let j = 0; j <= d; j++, x = mulMod(x, BigInt(k), mod)) {
                const valId = i * (d + 1) + j;
                M[k][valId] = x;
            }
        }
    }

    const { vec: C, sol } = systemOfLinearEquations(M, Array(N - r).fill(0n), mod);
    if (C.length === 0) return [];

    const res = Array.from({ length: r + 1 }, () => Array(d + 1).fill(0n));
    let noneZero = -2;
    {
        for (let i = 0; i <= d; i++) {
            if (sol[r * (d + 1) + i] !== 0n) {
                noneZero = -1;
                break;
            }
        }
        if (noneZero === -2) {
            for (let k = 0; k < C.length && noneZero === -2; k++) {
                for (let i = 0; i <= d; i++) {
                    if (C[k][r * (d + 1) + i] !== 0n) {
                        noneZero = k;
                        break;
                    }
                }
            }
        }
    }

    if (noneZero === -2) return [];

    for (let i = 0; i <= r; i++) {
        for (let j = 0; j <= d; j++) {
            res[i][j] = sol[i * (d + 1) + j];
            if (noneZero !== -1) {
                res[i][j] = addMod(res[i][j], C[noneZero][i * (d + 1) + j], mod);
            }
        }
    }
    return res;
}
// [0, R]
function random_int(R) {
    R++;
    return Math.floor(Math.random() * R);
}

// [L, R]
function random_int_range(L, R) {
    return random_int(R - L) + L;
}

// [L, R] * N
// O(N)
// O(NlogN) (is_unique)
function random_array(N, L, R, is_unique) {
    if (N == 0 || R < L) return [];
    if (is_unique) {
        if (R - L + 1 < N) return [];
        let res = [];
        let used = new Set();
        for (let i = 0; i < N; i++) {
            let x;
            while (true) {
                x = random_int(L, R);
                if (used.has(x)) continue;
                used.add(x);
                res.push(x);
                break;
            }
        }
        return res;
    } else {
        let res = [];
        for (let i = 0; i < N; i++) {
            res.push(random_int(L, R));
        }
        return res;
    }
}
// jsの数値型はdouble
// [0, R)
export function random_int(R) {
    return Math.floor(Math.random() * R);
}

// [0, R)
export function random_float(R) {
    return Math.random() * R;
}

// [L, R)
export function random_int_range(L, R) {
    let len = R - L;
    return random_int(len) + L;
}

// [L, R)
export function random_float_range(L, R) {
    let len = R - L;
    return random_float(len) + L;
}

// [L, R) * N
export function random_int_array(N, L, R) {
    let res = [];
    for (let i = 0; i < N; i++) {
        res.push(random_int_range(L, R));
    }
    return res;
}

// [L, R) * N
export function random_float_array(N, L, R) {
    let res = [];
    for (let i = 0; i < N; i++) {
        res.push(random_float_range(L, R));
    }
    return res;
}

export function random_shuffle(V) {
    let N = V.length, tmp, i;
    while (N) {
        i = Math.floor(Math.random() * N--);
        tmp = V[N];
        V[N] = V[i];
        V[i] = tmp;
    }
    return V;
}

// [1, N]
export function random_permutation(N) {
    let res = [];
    for (let i = 0; i < N; i++) {
        res.push(i + 1);
    }
    return random_shuffle(res);
}



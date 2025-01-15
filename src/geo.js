// python移植版でテスト
// https://judge.yosupo.jp/submission/261736

function cross(a, b, c) {
    return (b[1] - a[1]) * (c[0] - a[0]) - (c[1] - a[1]) * (b[0] - a[0]);
}

function lower_hull(P) {
    let N = P.length;
    let tmpP = [];
    for (let i = 0; i < N; i++) {
        tmpP.push([P[i][0], P[i][1], i]);
    }
    tmpP.sort((a, b) => { 
        if (a[0] == b[0]) return -a[1] + b[1];
        return a[0] - b[0]; 
    });
    let H = [];
    let len = 0;
    for (let i = 0; i < N; i++) {
        while (H.length >= 2) {
            len = H.length;
            if (cross(H[len - 2], H[len - 1], tmpP[i]) >= 0) {
                H.pop();
            } else {
                break;
            }
        }
        len = H.length;
        if (len == 0 || H[len - 1][0] != tmpP[i][0] || H[len - 1][1] != tmpP[i][1]) {
            H.push(tmpP[i]);
        }
    }
    len = H.length;
    if (len > 0 && H[0][0] != H[len - 1][0]) {
        for (let i = 0; i < N; i++) {
            if (H[len - 1][0] == tmpP[i][0] && H[len - 1][1] != tmpP[i][1]) {
                H.push(tmpP[i]);
                break;
            }
        }
    }
    return H;
}

function upper_hull(P) {
    let N = P.length;
    let tmpP = [];
    for (let i = 0; i < N; i++) {
        tmpP.push([P[i][0], P[i][1], i]);
    }
    tmpP.sort((a, b) => { 
        if (a[0] == b[0]) return a[1] - b[1];
        return a[0] - b[0]; 
    });
    let H = [];
    let len = 0;
    for (let i = 0; i < N; i++) {
        while (H.length >= 2) {
            len = H.length;
            if (cross(H[len - 2], H[len - 1], tmpP[i]) <= 0) {
                H.pop();
            } else {
                break;
            }
        }
        len = H.length;
        if (len == 0 || H[len - 1][0] != tmpP[i][0] || H[len - 1][1] != tmpP[i][1]) {
            H.push(tmpP[i]);
        }
    }
    len = H.length;
    if (len > 0 && H[0][0] != H[len - 1][0]) {
        for (let i = 0; i < N; i++) {
            if (H[len - 1][0] == tmpP[i][0] && H[len - 1][1] != tmpP[i][1]) {
                H.push(tmpP[i]);
                break;
            }
        }
    }
    return H;
}

function convex_hull(P) {
    let N = P.length;
    let tmpP = [];
    for (let i = 0; i < N; i++) {
        tmpP.push([P[i][0], P[i][1], i]);
    }
    tmpP.sort((a, b) => { 
        if (a[0] == b[0]) return -a[1] + b[1];
        return a[0] - b[0]; 
    });
    let len = 0;
    let HL = [], HU = [];
    for (let i = 0; i < N; i++) {
        while (HL.length >= 2) {
            len = HL.length;
            if (cross(HL[len - 2], HL[len - 1], tmpP[i]) >= 0) {
                HL.pop();
            } else {
                break;
            }
        }
        len = HL.length;
        if (len == 0 || HL[len - 1][0] != tmpP[i][0] || HL[len - 1][1] != tmpP[i][1]) {
            HL.push(tmpP[i]);
        }
    }
    for (let i = 0; i < N; i++) {
        while (HU.length >= 2) {
            len = HU.length;
            if (cross(HU[len - 2], HU[len - 1], tmpP[i]) <= 0) {
                HU.pop();
            } else {
                break;
            }
        }
        len = HU.length;
        if (len == 0 || HU[len - 1][0] != tmpP[i][0] || HU[len - 1][1] != tmpP[i][1]) {
            HU.push(tmpP[i]);
        }
    }

    len = HL.length;
    for (let i = HU.length - 1; i >= 0; i--) {
        if (len == 0 || (HL[0][0] != HU[i][0] || HL[0][1] != HU[i][1]) && (HL[len - 1][0] != HU[i][0] || HL[len - 1][1] != HU[i][1])) {
            HL.push(HU[i]);
        }
    }
    return HL
}
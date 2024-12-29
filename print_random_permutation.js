//import random_permutation from "./random.js"; できない

function random_shuffle(V) {
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
function random_permutation(N) {
    let res = [];
    for (let i = 0; i < N; i++) {
        res.push(i + 1);
    }
    return random_shuffle(res);
}

document.write("長さ10のランダムな順列を表示します。" + "<br>");

random_permutation(10).forEach(function(elem) {
    document.write(elem, " ");
});
document.write("<br>");